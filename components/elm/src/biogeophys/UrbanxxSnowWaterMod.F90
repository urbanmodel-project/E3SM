module UrbanxxSnowWaterMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Discrete multi-layer snowpack: Phase 2, SnowWater mass-balance and
  ! top-to-bottom meltwater percolation port. See
  ! plans/ELM_URBANxx_MultiLayerSnow_Port_Scoping_Plan.md.
  !
  ! Unlike UrbanxxSnowLayerMod's Phase 1 round-trip check (pure passthrough,
  ! no physics, seeded/read-back within the same call in SoilTemperatureMod),
  ! this module wires real physics with a NEW call site: ELM's own native
  ! SnowWater() runs in HydrologyNoDrainageMod.F90, AFTER Phasechange_beta
  ! has already modified snow-layer liquid/ice this same timestep (in
  ! SoilTemperatureMod.F90, where Phase 1's seed runs). So Phase 1's seed is
  ! stale by the time SnowWater needs input — this module re-seeds
  ! liquid/ice/dz fresh, immediately before ELM's own SnowWater() call, then
  ! calls URBANxx's own SnowWater port on that freshly-seeded state. The
  ! check (called immediately after ELM's native SnowWater()) compares
  ! URBANxx's evolved state against ELM's own newly-evolved state — a real
  ! physics-parity comparison to a roundoff tolerance, not a bit-identity
  ! check.
  !
  ! NumberOfActiveSnowLayers is NOT re-seeded here: snl(c) does not change
  ! between SoilTemperatureMod (Phase 1's seed point) and this module's call
  ! site (layer-count changes only happen later, via CombineSnowLayers/
  ! DivideSnowLayers, which run after SnowWater — Phase 3), so Phase 1's
  ! already-seeded value is still correct.

  use iso_c_binding
  use urban_mod
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use spmdMod            , only : masterproc, iam
  use elm_varctl         , only : iulog
  use abortutils         , only : endrun
  use UrbanxxInstanceMod , only : urbanxx
  use UrbanxxSnowLayerMod, only : PackLayerQuantity, max_snow_layers

  implicit none

  private

  public :: urbanxx_snowWater
  public :: urbanxx_snowWater_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowWater(num_urbanl, num_urbanc, filter_urbanc, begc, dtime)
    !
    ! !DESCRIPTION:
    ! Re-seed URBANxx's snow-layer liquid water, ice content, and thickness
    ! from ELM's current (post-Phasechange_beta, pre-native-SnowWater)
    ! state, then run URBANxx's own SnowWater mass-balance/percolation.
    !
    use ColumnDataType         , only : col_ws
    use ColumnType             , only : col_pp
    use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer        , intent(in) :: filter_urbanc(:)
    integer        , intent(in) :: begc
    real(r8)       , intent(in) :: dtime

    logical :: isLayoutLeft
    real(c_double), allocatable, target :: buf_roof(:), buf_improad(:), buf_perroad(:)
    integer(c_int) :: status
    integer(c_int), dimension(2) :: size2D

    isLayoutLeft = UrbanKokkosIsLayoutLeft()
    size2D(1) = num_urbanl
    size2D(2) = max_snow_layers

    allocate(buf_roof(num_urbanl * max_snow_layers))
    allocate(buf_improad(num_urbanl * max_snow_layers))
    allocate(buf_perroad(num_urbanl * max_snow_layers))

    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_ws%h2osoi_liq, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call UrbanSetSnowLayerLiquidWaterRoof(urbanxx, c_loc(buf_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerLiquidWaterImperviousRoad(urbanxx, c_loc(buf_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerLiquidWaterPerviousRoad(urbanxx, c_loc(buf_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_ws%h2osoi_ice, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call UrbanSetSnowLayerIceContentRoof(urbanxx, c_loc(buf_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerIceContentImperviousRoad(urbanxx, c_loc(buf_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerIceContentPerviousRoad(urbanxx, c_loc(buf_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_pp%dz, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call UrbanSetSnowLayerThicknessRoof(urbanxx, c_loc(buf_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerThicknessImperviousRoad(urbanxx, c_loc(buf_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerThicknessPerviousRoad(urbanxx, c_loc(buf_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    deallocate(buf_roof, buf_improad, buf_perroad)

    call UrbanComputeSnowWater(urbanxx, real(dtime, c_double), status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine urbanxx_snowWater

  !-----------------------------------------------------------------------
  subroutine CheckOneLayerQuantity(label, urbanxx_buf, elm_buf, n, tol)
    !
    ! !DESCRIPTION:
    ! Compare URBANxx's evolved snow-layer state against ELM's own
    ! independently-evolved state (both post their respective SnowWater
    ! calls this timestep). Real floating-point physics on both sides, so
    ! compared to a roundoff tolerance rather than bit-identity.
    !
    implicit none
    character(len=*), intent(in) :: label
    real(c_double)  , intent(in) :: urbanxx_buf(:)
    real(r8)        , intent(in) :: elm_buf(:)
    integer          , intent(in) :: n
    real(r8)        , intent(in) :: tol

    integer :: i
    real(r8) :: max_err

    max_err = 0._r8
    do i = 1, n
      max_err = max(max_err, abs(real(urbanxx_buf(i), r8) - elm_buf(i)))
    end do
    if (max_err > tol) then
      write(iulog,'(A,A,A,ES24.16)') 'ERROR: UrbanxxSnowWaterMod mismatch for ', &
           trim(label), ', max_err = ', max_err
      call endrun(msg='UrbanxxSnowWaterMod: SnowWater physics-parity check failed for '//trim(label))
    else
      write(iulog,'(A,A,A,ES24.16)') 'Max error in UrbanxxSnowWaterMod ', trim(label), ' : ', max_err
    end if

  end subroutine CheckOneLayerQuantity

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowWater_check(num_urbanl, num_urbanc, filter_urbanc, begc)
    !
    ! !DESCRIPTION:
    ! Called immediately after ELM's native SnowWater() in
    ! HydrologyNoDrainageMod.F90. Reads back URBANxx's evolved
    ! liquid/ice/dz and compares against ELM's freshly-updated values.
    !
    use ColumnDataType         , only : col_ws
    use ColumnType             , only : col_pp
    use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer        , intent(in) :: filter_urbanc(:)
    integer        , intent(in) :: begc

    logical :: isLayoutLeft
    real(c_double), allocatable, target :: urbanxx_roof(:), urbanxx_improad(:), urbanxx_perroad(:)
    real(c_double), allocatable, target :: elm_roof(:), elm_improad(:), elm_perroad(:)
    integer(c_int) :: status
    integer(c_int), dimension(2) :: size2D
    integer :: n

    isLayoutLeft = UrbanKokkosIsLayoutLeft()
    size2D(1) = num_urbanl
    size2D(2) = max_snow_layers
    n = num_urbanl * max_snow_layers

    allocate(urbanxx_roof(n)); allocate(urbanxx_improad(n)); allocate(urbanxx_perroad(n))
    allocate(elm_roof(n)); allocate(elm_improad(n)); allocate(elm_perroad(n))

    ! Liquid water: tolerance matches the 1e-9 standard used by other
    ! Urbanxx*_check routines (e.g. urbanxx_soilTemperature_check).
    call UrbanGetSnowLayerLiquidWaterRoof(urbanxx, c_loc(urbanxx_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerLiquidWaterImperviousRoad(urbanxx, c_loc(urbanxx_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerLiquidWaterPerviousRoad(urbanxx, c_loc(urbanxx_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_ws%h2osoi_liq, isLayoutLeft, elm_roof, elm_improad, elm_perroad)
    call CheckOneLayerQuantity('SnowLayerLiquidWaterRoof', urbanxx_roof, elm_roof, n, 1.0e-9_r8)
    call CheckOneLayerQuantity('SnowLayerLiquidWaterImperviousRoad', urbanxx_improad, elm_improad, n, 1.0e-9_r8)
    call CheckOneLayerQuantity('SnowLayerLiquidWaterPerviousRoad', urbanxx_perroad, elm_perroad, n, 1.0e-9_r8)

    call UrbanGetSnowLayerIceContentRoof(urbanxx, c_loc(urbanxx_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerIceContentImperviousRoad(urbanxx, c_loc(urbanxx_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerIceContentPerviousRoad(urbanxx, c_loc(urbanxx_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_ws%h2osoi_ice, isLayoutLeft, elm_roof, elm_improad, elm_perroad)
    call CheckOneLayerQuantity('SnowLayerIceContentRoof', urbanxx_roof, elm_roof, n, 1.0e-9_r8)
    call CheckOneLayerQuantity('SnowLayerIceContentImperviousRoad', urbanxx_improad, elm_improad, n, 1.0e-9_r8)
    call CheckOneLayerQuantity('SnowLayerIceContentPerviousRoad', urbanxx_perroad, elm_perroad, n, 1.0e-9_r8)

    ! Thickness: much smaller magnitude (meters), tighter absolute tolerance.
    call UrbanGetSnowLayerThicknessRoof(urbanxx, c_loc(urbanxx_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerThicknessImperviousRoad(urbanxx, c_loc(urbanxx_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerThicknessPerviousRoad(urbanxx, c_loc(urbanxx_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_pp%dz, isLayoutLeft, elm_roof, elm_improad, elm_perroad)
    call CheckOneLayerQuantity('SnowLayerThicknessRoof', urbanxx_roof, elm_roof, n, 1.0e-12_r8)
    call CheckOneLayerQuantity('SnowLayerThicknessImperviousRoad', urbanxx_improad, elm_improad, n, 1.0e-12_r8)
    call CheckOneLayerQuantity('SnowLayerThicknessPerviousRoad', urbanxx_perroad, elm_perroad, n, 1.0e-12_r8)

    deallocate(urbanxx_roof, urbanxx_improad, urbanxx_perroad)
    deallocate(elm_roof, elm_improad, elm_perroad)

  end subroutine urbanxx_snowWater_check

end module UrbanxxSnowWaterMod
