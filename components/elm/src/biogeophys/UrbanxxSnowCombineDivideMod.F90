module UrbanxxSnowCombineDivideMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Discrete multi-layer snowpack: Phase 3, CombineSnowLayers/DivideSnowLayers
  ! port (layer-count evolution: merging thin layers, splitting thick ones).
  ! See plans/ELM_URBANxx_MultiLayerSnow_Port_Scoping_Plan.md.
  !
  ! Like UrbanxxSnowCompactionMod (Phase 2), this does NOT re-seed
  ! temperature/liquid/ice/dz/NumberOfActiveSnowLayers: URBANxx's own
  ! SnowWater + SnowCompaction ports already evolved that state earlier this
  ! same timestep, and nothing between SnowCompaction and this call site
  ! modifies it. Both ELM's CombineSnowLayers and DivideSnowLayers operate
  ! purely on state already owned by whichever side is running them (no
  ! external forcing, no soil-side or atmosphere-side inputs) — merges
  ! average two existing layers' state together, and splits derive a new
  ! layer's state entirely from the parent layer being split. This resolves
  ! the scoping doc's "Open sub-question for Phase 1": no event-driven
  ! reseed from ELM is needed for Combine/Divide-driven layer-count changes.
  ! (A brand-new layer created by fresh snowfall, i.e. ELM's snl 0 -> -1
  ! transition in CanopyHydrologyMod.F90, is a separate call site entirely
  ! and is out of scope here — see the scoping doc's "Important scoping
  ! correction": that event is what the still-open Known Gap #1 abort is
  ! about, and closing it is Phase 4's job, not this one.)
  !
  ! NumberOfActiveSnowLayers becomes an OUTPUT of this module for the first
  ! time in the port (Phases 1-2 only ever read/re-seeded it) — URBANxx's
  ! own Combine/Divide port computes the new active-layer count internally.

  use iso_c_binding
  use urban_mod
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use spmdMod            , only : masterproc, iam
  use elm_varctl         , only : iulog
  use abortutils         , only : endrun
  use UrbanxxInstanceMod , only : urbanxx
  use UrbanxxSnowLayerMod, only : PackLayerQuantity, PackActiveCount, max_snow_layers
  use UrbanxxDiagnosticControlMod, only : urbanxx_softfail_permits, urbanxx_packed_tol, &
       TOL_PHYS_TEMP, TOL_PHYS_SNOW_WATER, TOL_PHYS_SNOW_DZ

  implicit none

  private

  public :: urbanxx_snowCombineDivide
  public :: urbanxx_snowCombineDivide_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowCombineDivide(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Run URBANxx's own CombineSnowLayers, then DivideSnowLayers, continuing
    ! from the liquid/ice/dz/active-count state ComputeSnowCompaction already
    ! evolved earlier this timestep. No inputs need seeding.
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl

    integer(c_int) :: status

    call UrbanComputeSnowLayerCombine(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanComputeSnowLayerDivide(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Phase 6 fix (netLongwave-mismatch-at-timestep-124 root cause): neither
    ! SnowCompaction nor Combine/Divide above ever touch the scalar
    ! SnowDepth view URBANxx's ComputeSnowCover/ComputeUpdateSnowFraction
    ! read for FracSno — resync it here from the just-restructured layer
    ! thicknesses, mirroring ELM's own snow_depth(c) recompute inside
    ! CombineSnowLayers. See UrbanSnowCombineLayers.cpp for the full
    ! explanation.
    call UrbanSyncSnowDepthFromLayers(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine urbanxx_snowCombineDivide

  !-----------------------------------------------------------------------
  subroutine CheckActiveCount(label, urbanxx_buf, elm_buf, n)
    !
    ! !DESCRIPTION:
    ! Exact (integer) comparison of the new active-layer count.
    !
    implicit none
    character(len=*), intent(in) :: label
    real(c_double)  , intent(in) :: urbanxx_buf(:)
    real(c_double)  , intent(in) :: elm_buf(:)
    integer          , intent(in) :: n

    integer :: i
    real(r8) :: max_err

    max_err = 0._r8
    do i = 1, n
      max_err = max(max_err, abs(real(urbanxx_buf(i), r8) - real(elm_buf(i), r8)))
    end do
    if (max_err > 0._r8) then
      write(iulog,'(A,A,A,ES24.16)') 'ERROR: UrbanxxSnowCombineDivideMod mismatch for ', &
           trim(label), ', max_err = ', max_err
      if (urbanxx_softfail_permits(max_err, 0._r8)) then
        write(iulog,*) '  URBANxx check failure permitted to continue: '//trim(label)// &
             ' active-layer-count mismatch (urbanxx_softfail/urbanxx_always_continue); continuing'
      else
        call endrun(msg='UrbanxxSnowCombineDivideMod: active-layer-count parity check failed for '//trim(label))
      end if
    else
      write(iulog,'(A,A,A,ES24.16)') 'Max error in UrbanxxSnowCombineDivideMod ', trim(label), ' : ', max_err
    end if

  end subroutine CheckActiveCount

  !-----------------------------------------------------------------------
  subroutine CheckOneLayerQuantity(label, urbanxx_buf, elm_buf, n, tol, tolv)
    implicit none
    character(len=*), intent(in) :: label
    real(c_double)  , intent(in) :: urbanxx_buf(:)
    real(r8)        , intent(in) :: elm_buf(:)
    integer          , intent(in) :: n
    real(r8)        , intent(in) :: tol
    real(r8)        , intent(in) :: tolv(:)   ! per-element tolerance (two-tier)

    integer :: i
    real(r8) :: max_err, max_ratio

    max_err = 0._r8
    max_ratio = 0._r8
    do i = 1, n
      max_err = max(max_err, abs(real(urbanxx_buf(i), r8) - elm_buf(i)))
      max_ratio = max(max_ratio, abs(real(urbanxx_buf(i), r8) - elm_buf(i)) / tolv(i))
    end do
    if (max_ratio > 1._r8) then
      write(iulog,'(A,A,A,ES24.16)') 'ERROR: UrbanxxSnowCombineDivideMod mismatch for ', &
           trim(label), ', max_err = ', max_err
      if (urbanxx_softfail_permits(max_ratio, 1._r8)) then
        write(iulog,*) '  URBANxx check failure permitted to continue (urbanxx_softfail/urbanxx_always_continue); logging and continuing'
      else
        call endrun(msg='UrbanxxSnowCombineDivideMod: Combine/Divide physics-parity check failed for '//trim(label))
      end if
    else
      write(iulog,'(A,A,A,ES24.16)') 'Max error in UrbanxxSnowCombineDivideMod ', trim(label), ' : ', max_err
    end if

  end subroutine CheckOneLayerQuantity

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowCombineDivide_check(num_urbanl, num_urbanc, filter_urbanc, begc)
    !
    ! !DESCRIPTION:
    ! Called immediately after ELM's native CombineSnowLayers()/
    ! DivideSnowLayers() in HydrologyNoDrainageMod.F90 (before the "set
    ! empty snow layers to zero" step — URBANxx's own Divide port already
    ! zero-fills inactive slots itself). Compares URBANxx's new active-layer
    ! count and evolved layer state against ELM's freshly restructured
    ! state.
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
    real(r8), allocatable :: tolT(:), tolW(:), tolZ(:)

    isLayoutLeft = UrbanKokkosIsLayoutLeft()
    size2D(1) = num_urbanl
    size2D(2) = max_snow_layers

    ! --- Active-layer count (exact match expected) ---
    allocate(urbanxx_roof(num_urbanl)); allocate(urbanxx_improad(num_urbanl)); allocate(urbanxx_perroad(num_urbanl))
    allocate(elm_roof(num_urbanl)); allocate(elm_improad(num_urbanl)); allocate(elm_perroad(num_urbanl))

    call UrbanGetNumberOfActiveSnowLayersRoof(urbanxx, c_loc(urbanxx_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetNumberOfActiveSnowLayersImperviousRoad(urbanxx, c_loc(urbanxx_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetNumberOfActiveSnowLayersPerviousRoad(urbanxx, c_loc(urbanxx_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call PackActiveCount(num_urbanc, filter_urbanc, elm_roof, elm_improad, elm_perroad)
    call CheckActiveCount('NumberOfActiveSnowLayersRoof', urbanxx_roof, elm_roof, num_urbanl)
    call CheckActiveCount('NumberOfActiveSnowLayersImperviousRoad', urbanxx_improad, elm_improad, num_urbanl)
    call CheckActiveCount('NumberOfActiveSnowLayersPerviousRoad', urbanxx_perroad, elm_perroad, num_urbanl)

    deallocate(urbanxx_roof, urbanxx_improad, urbanxx_perroad)
    deallocate(elm_roof, elm_improad, elm_perroad)

    ! --- Layer temperature, liquid water, ice content, thickness ---
    n = num_urbanl * max_snow_layers
    allocate(tolT(n), tolW(n), tolZ(n))
    call urbanxx_packed_tol(num_urbanc, filter_urbanc, num_urbanl, max_snow_layers, &
         isLayoutLeft, 1.0e-9_r8, TOL_PHYS_TEMP, tolT)
    call urbanxx_packed_tol(num_urbanc, filter_urbanc, num_urbanl, max_snow_layers, &
         isLayoutLeft, 1.0e-9_r8, TOL_PHYS_SNOW_WATER, tolW)
    call urbanxx_packed_tol(num_urbanc, filter_urbanc, num_urbanl, max_snow_layers, &
         isLayoutLeft, 1.0e-12_r8, TOL_PHYS_SNOW_DZ, tolZ)
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
    call CheckOneLayerQuantity('SnowLayerTempRoof', urbanxx_roof, elm_roof, n, 1.0e-9_r8, tolT)
    call CheckOneLayerQuantity('SnowLayerTempImperviousRoad', urbanxx_improad, elm_improad, n, 1.0e-9_r8, tolT)
    call CheckOneLayerQuantity('SnowLayerTempPerviousRoad', urbanxx_perroad, elm_perroad, n, 1.0e-9_r8, tolT)

    call UrbanGetSnowLayerLiquidWaterRoof(urbanxx, c_loc(urbanxx_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerLiquidWaterImperviousRoad(urbanxx, c_loc(urbanxx_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerLiquidWaterPerviousRoad(urbanxx, c_loc(urbanxx_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_ws%h2osoi_liq, isLayoutLeft, elm_roof, elm_improad, elm_perroad)
    call CheckOneLayerQuantity('SnowLayerLiquidWaterRoof', urbanxx_roof, elm_roof, n, 1.0e-9_r8, tolW)
    call CheckOneLayerQuantity('SnowLayerLiquidWaterImperviousRoad', urbanxx_improad, elm_improad, n, 1.0e-9_r8, tolW)
    call CheckOneLayerQuantity('SnowLayerLiquidWaterPerviousRoad', urbanxx_perroad, elm_perroad, n, 1.0e-9_r8, tolW)

    call UrbanGetSnowLayerIceContentRoof(urbanxx, c_loc(urbanxx_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerIceContentImperviousRoad(urbanxx, c_loc(urbanxx_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerIceContentPerviousRoad(urbanxx, c_loc(urbanxx_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_ws%h2osoi_ice, isLayoutLeft, elm_roof, elm_improad, elm_perroad)
    call CheckOneLayerQuantity('SnowLayerIceContentRoof', urbanxx_roof, elm_roof, n, 1.0e-9_r8, tolW)
    call CheckOneLayerQuantity('SnowLayerIceContentImperviousRoad', urbanxx_improad, elm_improad, n, 1.0e-9_r8, tolW)
    call CheckOneLayerQuantity('SnowLayerIceContentPerviousRoad', urbanxx_perroad, elm_perroad, n, 1.0e-9_r8, tolW)

    call UrbanGetSnowLayerThicknessRoof(urbanxx, c_loc(urbanxx_roof), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerThicknessImperviousRoad(urbanxx, c_loc(urbanxx_improad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSnowLayerThicknessPerviousRoad(urbanxx, c_loc(urbanxx_perroad), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_pp%dz, isLayoutLeft, elm_roof, elm_improad, elm_perroad)
    call CheckOneLayerQuantity('SnowLayerThicknessRoof', urbanxx_roof, elm_roof, n, 1.0e-12_r8, tolZ)
    call CheckOneLayerQuantity('SnowLayerThicknessImperviousRoad', urbanxx_improad, elm_improad, n, 1.0e-12_r8, tolZ)
    call CheckOneLayerQuantity('SnowLayerThicknessPerviousRoad', urbanxx_perroad, elm_perroad, n, 1.0e-12_r8, tolZ)

    deallocate(urbanxx_roof, urbanxx_improad, urbanxx_perroad)
    deallocate(elm_roof, elm_improad, elm_perroad)

  end subroutine urbanxx_snowCombineDivide_check

end module UrbanxxSnowCombineDivideMod
