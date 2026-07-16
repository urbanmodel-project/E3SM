module UrbanxxSnowCompactionMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Discrete multi-layer snowpack: Phase 2, SnowCompaction port (layer
  ! thickness evolution due to destructive metamorphism, overburden, and
  ! melt). See plans/ELM_URBANxx_MultiLayerSnow_Port_Scoping_Plan.md.
  !
  ! Unlike urbanxx_snowWater (UrbanxxSnowWaterMod.F90), this does NOT
  ! re-seed liquid/ice/dz: URBANxx's own SnowWater port already evolved
  ! those views earlier this same timestep (called from the SnowWater call
  ! site in HydrologyNoDrainageMod.F90), and nothing between that call site
  ! and this one (SurfaceRunoff/Infiltration/SoilWater/WaterTable) modifies
  ! snow-layer liquid/ice/dz — only soil water/drainage state. So this
  ! module continues directly from where ComputeSnowWater left off.
  !
  ! Two NEW per-layer inputs are seeded here: imelt (col_ef%imelt, set by
  ! Phasechange_beta in SoilTemperatureMod.F90, step 3 — already finalized
  ! for this timestep by the time this runs in step 5) and frac_iceold
  ! (col_ws%frac_iceold, snapshotted at the very start of the timestep in
  ! elm_drv_init, before any physics this timestep — also stable and valid
  ! by this point).

  use iso_c_binding
  use urban_mod
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use spmdMod            , only : masterproc, iam
  use elm_varctl         , only : iulog
  use abortutils         , only : endrun
  use UrbanxxInstanceMod , only : urbanxx
  use UrbanxxSnowLayerMod, only : PackLayerQuantity, PackLayerQuantityInt, max_snow_layers
  use UrbanxxDiagnosticControlMod, only : urbanxx_softfail_permits

  implicit none

  private

  public :: urbanxx_snowCompaction
  public :: urbanxx_snowCompaction_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowCompaction(num_urbanl, num_urbanc, filter_urbanc, begc, dtime)
    !
    ! !DESCRIPTION:
    ! Seed imelt and frac_iceold, then run URBANxx's own SnowCompaction on
    ! the snow-layer state already evolved by ComputeSnowWater this
    ! timestep.
    !
    use ColumnDataType         , only : col_ws, col_ef
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

    call PackLayerQuantityInt(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_ef%imelt, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call UrbanSetSnowLayerImeltRoof(urbanxx, c_loc(buf_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerImeltImperviousRoad(urbanxx, c_loc(buf_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerImeltPerviousRoad(urbanxx, c_loc(buf_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_ws%frac_iceold, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call UrbanSetSnowLayerFracIceOldRoof(urbanxx, c_loc(buf_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerFracIceOldImperviousRoad(urbanxx, c_loc(buf_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetSnowLayerFracIceOldPerviousRoad(urbanxx, c_loc(buf_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    deallocate(buf_roof, buf_improad, buf_perroad)

    call UrbanComputeSnowCompaction(urbanxx, real(dtime, c_double), status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine urbanxx_snowCompaction

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
      write(iulog,'(A,A,A,ES24.16)') 'ERROR: UrbanxxSnowCompactionMod mismatch for ', &
           trim(label), ', max_err = ', max_err
      if (urbanxx_softfail_permits(max_err, tol)) then
        write(iulog,*) '  URBANxx check failure permitted to continue (urbanxx_softfail/urbanxx_always_continue); logging and continuing'
      else
        call endrun(msg='UrbanxxSnowCompactionMod: SnowCompaction physics-parity check failed for '//trim(label))
      end if
    else
      write(iulog,'(A,A,A,ES24.16)') 'Max error in UrbanxxSnowCompactionMod ', trim(label), ' : ', max_err
    end if

  end subroutine CheckOneLayerQuantity

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowCompaction_check(num_urbanl, num_urbanc, filter_urbanc, begc)
    !
    ! !DESCRIPTION:
    ! Called immediately after ELM's native SnowCompaction() in
    ! HydrologyNoDrainageMod.F90. Compares URBANxx's evolved dz against
    ! ELM's freshly-updated dz.
    !
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

  end subroutine urbanxx_snowCompaction_check

end module UrbanxxSnowCompactionMod
