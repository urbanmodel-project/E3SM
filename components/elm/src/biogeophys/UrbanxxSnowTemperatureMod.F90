module UrbanxxSnowTemperatureMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Discrete multi-layer snowpack: Phase 4, multi-layer temperature solve
  ! seed + check. See plans/ELM_URBANxx_MultiLayerSnow_Port_Scoping_Plan.md.
  !
  ! URBANxx's snow-layer temperature, liquid water, and ice content are NOT
  ! reseeded from ELM every timestep (unlike Phase 1's original design for
  ! these same quantities) — once ComputeHeatDiffusion (via
  ! urbanxx_soilTemperature) actually evolves them, re-seeding afterward
  ! would silently discard the just-solved result. Instead:
  !
  !   1. urbanxx_snowTemperatureSeed seeds NumberOfActiveSnowLayers, Temp,
  !      Liq, Ice, and Dz fresh from ELM's CURRENT state, called BEFORE
  !      urbanxx_soilTemperature (a new call site in SoilTemperatureMod.F90,
  !      earlier than Phase 1's UrbanxxSnowLayerMod%urbanxx_snowLayer, which
  !      still separately re-seeds NumberOfActiveSnowLayers/Dz/Zc afterward
  !      for its own pure-plumbing round-trip check — harmless duplication,
  !      since layer count and thickness don't change during the solve).
  !      This is what lets ComputeHeatDiffusion's combined snow+soil solve
  !      see a snow layer CanopyHydrology created earlier this same
  !      timestep, rather than one timestep stale.
  !   2. ComputeHeatDiffusion (Phase 4c/4d) evolves SnowLayerTemp and
  !      phase-change-adjusts SnowLayerIceContent/LiquidWater.
  !   3. urbanxx_snowTemperature_check compares URBANxx's post-solve values
  !      against ELM's independently post-solve values (SolveTemperature +
  !      Phasechange_beta) — a real physics-parity check to a roundoff
  !      tolerance, analogous to Phase 2/3's checks, not a bit-identical
  !      round-trip. Called after Phasechange_beta in SoilTemperatureMod.F90,
  !      alongside the existing urbanxx_soilTemperature_check.

  use iso_c_binding
  use urban_mod
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use spmdMod            , only : masterproc, iam
  use elm_varctl         , only : iulog
  use abortutils         , only : endrun
  use UrbanxxInstanceMod , only : urbanxx
  use UrbanxxSnowLayerMod, only : PackLayerQuantity, PackActiveCount, max_snow_layers

  implicit none

  private

  public :: urbanxx_snowTemperatureSeed
  public :: urbanxx_snowTemperature_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowTemperatureSeed(num_urbanl, num_urbanc, filter_urbanc, begc)
    !
    ! !DESCRIPTION:
    ! Seed NumberOfActiveSnowLayers, layer temperature, liquid water, ice
    ! content, and thickness from ELM's current state, immediately before
    ! urbanxx_soilTemperature (and its internal ComputeHeatDiffusion) runs.
    !
    use ColumnDataType         , only : col_es, col_ws
    use ColumnType             , only : col_pp
    use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer        , intent(in) :: filter_urbanc(:)
    integer        , intent(in) :: begc

    logical :: isLayoutLeft
    real(c_double), allocatable, target :: buf_roof(:), buf_improad(:), buf_perroad(:)
    integer(c_int) :: status
    integer(c_int), dimension(2) :: size2D

    isLayoutLeft = UrbanKokkosIsLayoutLeft()
    size2D(1) = num_urbanl
    size2D(2) = max_snow_layers

    allocate(buf_roof(num_urbanl))
    allocate(buf_improad(num_urbanl))
    allocate(buf_perroad(num_urbanl))
    call PackActiveCount(num_urbanc, filter_urbanc, buf_roof, buf_improad, buf_perroad)
    call UrbanSetNumberOfActiveSnowLayersRoof(urbanxx, c_loc(buf_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetNumberOfActiveSnowLayersImperviousRoad(urbanxx, c_loc(buf_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetNumberOfActiveSnowLayersPerviousRoad(urbanxx, c_loc(buf_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    deallocate(buf_roof, buf_improad, buf_perroad)

    allocate(buf_roof(num_urbanl * max_snow_layers))
    allocate(buf_improad(num_urbanl * max_snow_layers))
    allocate(buf_perroad(num_urbanl * max_snow_layers))

    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_es%t_soisno, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call UrbanSetSnowLayerTempRoof(urbanxx, c_loc(buf_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerTempImperviousRoad(urbanxx, c_loc(buf_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerTempPerviousRoad(urbanxx, c_loc(buf_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

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

  end subroutine urbanxx_snowTemperatureSeed

  !-----------------------------------------------------------------------
  subroutine CheckOneLayerQuantity(label, urbanxx_buf, elm_buf, n, tol)
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
      write(iulog,'(A,A,A,ES24.16)') 'ERROR: UrbanxxSnowTemperatureMod mismatch for ', &
           trim(label), ', max_err = ', max_err
      call endrun(msg='UrbanxxSnowTemperatureMod: multi-layer temperature-solve parity check failed for '//trim(label))
    else
      write(iulog,'(A,A,A,ES24.16)') 'Max error in UrbanxxSnowTemperatureMod ', trim(label), ' : ', max_err
    end if

  end subroutine CheckOneLayerQuantity

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowTemperature_check(num_urbanl, num_urbanc, filter_urbanc, begc)
    !
    ! !DESCRIPTION:
    ! Called after ELM's native SolveTemperature + Phasechange_beta in
    ! SoilTemperatureMod.F90. Compares URBANxx's post-solve snow-layer
    ! temperature/liquid/ice against ELM's freshly-updated values.
    !
    use ColumnDataType         , only : col_es, col_ws
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

    call UrbanGetSnowLayerTempRoof(urbanxx, c_loc(urbanxx_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerTempImperviousRoad(urbanxx, c_loc(urbanxx_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerTempPerviousRoad(urbanxx, c_loc(urbanxx_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_es%t_soisno, isLayoutLeft, elm_roof, elm_improad, elm_perroad)
    call CheckOneLayerQuantity('SnowLayerTempRoof', urbanxx_roof, elm_roof, n, 1.0e-9_r8)
    call CheckOneLayerQuantity('SnowLayerTempImperviousRoad', urbanxx_improad, elm_improad, n, 1.0e-9_r8)
    call CheckOneLayerQuantity('SnowLayerTempPerviousRoad', urbanxx_perroad, elm_perroad, n, 1.0e-9_r8)

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

    deallocate(urbanxx_roof, urbanxx_improad, urbanxx_perroad)
    deallocate(elm_roof, elm_improad, elm_perroad)

  end subroutine urbanxx_snowTemperature_check

end module UrbanxxSnowTemperatureMod
