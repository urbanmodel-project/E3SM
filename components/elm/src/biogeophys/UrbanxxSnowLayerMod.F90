module UrbanxxSnowLayerMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Discrete multi-layer snowpack: Phase 1 data model + plumbing (see
  ! plans/ELM_URBANxx_MultiLayerSnow_Port_Scoping_Plan.md), trimmed down by
  ! Phase 4 (multi-layer temperature solve).
  !
  ! This module now wires only 3 of the 9 planned quantities as a pure,
  ! every-timestep round-trip check with no physics in between:
  ! NumberOfActiveSnowLayers (unaffected by the Phase 4 solve — layer-count
  ! changes only happen in HydrologyNoDrainageMod.F90's Combine/Divide,
  ! Phase 3, which runs later in the timestep), Dz (also unaffected — phase
  ! change modifies temperature/liquid/ice, never thickness), and Zc (kept
  ! for backward-compatible plumbing validation, though the Phase 4 solve
  ! computes its own transient snow Zc/Zi internally rather than reading
  ! this view — see SolveHeatDiffusionSnowSoil in UrbanHeatDiffusion.cpp).
  !
  ! Temperature, liquid water, and ice content moved OUT of this module's
  ! round-trip seed+check: once Phase 4's ComputeHeatDiffusion actually
  ! evolves SnowLayerTemp (and phase-change-adjusts liquid/ice) instead of
  ! leaving those views untouched, seeding them here (which used to run
  ! AFTER urbanxx_soilTemperature) would silently overwrite URBANxx's
  ! just-solved values with ELM's pre-solve state every timestep, defeating
  ! the whole point of the solve. See UrbanxxSnowTemperatureMod.F90, which
  ! now (a) seeds those 3 quantities fresh from ELM BEFORE
  ! urbanxx_soilTemperature runs (so the solve gets current input, including
  ! any snow layer CanopyHydrology created earlier this same timestep), and
  ! (b) checks URBANxx's post-solve values against ELM's independently
  ! post-solve values (a real physics-parity check, not a bit-identical
  ! round-trip).
  !
  ! Layer-interface depth (Zi) and thermal conductivity/heat capacity remain
  ! NOT wired here (Phase 4 confirmed they should never be — ELM's own tk/cv
  ! are transient locals never exposed outside SoilTemperature's call chain,
  ! so URBANxx computes its own transiently too, entirely from state it
  ! already owns; there is nothing to seed).
  !
  ! Indexing convention: URBANxx cpp_index 0 = snow layer nearest the soil,
  ! increasing away from soil. Maps to ELM's snl+1..0 range via a plain sign
  ! flip: elm_index = -(cpp_index). See UrbanSnowLayerSetters.cpp.

  use iso_c_binding
  use urban_mod
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use spmdMod            , only : masterproc, iam
  use elm_varctl         , only : iulog
  use abortutils         , only : endrun
  use UrbanxxInstanceMod , only : urbanxx

  implicit none

  private

  ! Must match URBANXX::NUM_SNOW_LAYERS (UrbanConstants.h). Checked at
  ! runtime against ELM's own nlevsno in urbanxx_snowLayer_init.
  integer, parameter :: max_snow_layers = 5

  ! Persistent round-trip (readback) buffers, allocated once in init.
  ! "_rb" = value read back from URBANxx immediately after being set —
  ! compared against the freshly-packed ELM value in urbanxx_snowLayer_check.
  real(c_double), allocatable, target :: numactive_roof_rb(:)
  real(c_double), allocatable, target :: numactive_improad_rb(:)
  real(c_double), allocatable, target :: numactive_perroad_rb(:)

  real(c_double), allocatable, target :: dz_roof_rb(:)
  real(c_double), allocatable, target :: dz_improad_rb(:)
  real(c_double), allocatable, target :: dz_perroad_rb(:)

  real(c_double), allocatable, target :: zc_roof_rb(:)
  real(c_double), allocatable, target :: zc_improad_rb(:)
  real(c_double), allocatable, target :: zc_perroad_rb(:)

  public :: urbanxx_snowLayer_init
  public :: urbanxx_snowLayer
  public :: urbanxx_snowLayer_check
  ! Exposed for reuse by UrbanxxSnowWaterMod.F90 (Phase 2), which needs to
  ! re-seed liquid/ice/dz at a later point in the timestep (after
  ! Phasechange_beta, right before ELM's native SnowWater()) than this
  ! module's own every-timestep Phase 1 seed (which runs earlier, in
  ! SoilTemperatureMod.F90, and is therefore stale by the time SnowWater
  ! needs fresh input).
  public :: PackLayerQuantity
  public :: PackLayerQuantityInt
  public :: PackActiveCount
  public :: max_snow_layers

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowLayer_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent round-trip buffers. Called once during
    ! initialization.
    !
    use elm_varpar, only : nlevsno
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl

    if (nlevsno /= max_snow_layers) then
      call endrun(msg='UrbanxxSnowLayerMod: nlevsno /= max_snow_layers (5). '// &
           'Phase 1 of the multi-layer snow port assumes use_extrasnowlayers=.false. '// &
           '(nlevsno=5); see ELM_URBANxx_MultiLayerSnow_Port_Scoping_Plan.md Open Questions.')
    end if

    allocate(numactive_roof_rb(num_urbanl))
    allocate(numactive_improad_rb(num_urbanl))
    allocate(numactive_perroad_rb(num_urbanl))

    allocate(dz_roof_rb(num_urbanl * max_snow_layers))
    allocate(dz_improad_rb(num_urbanl * max_snow_layers))
    allocate(dz_perroad_rb(num_urbanl * max_snow_layers))

    allocate(zc_roof_rb(num_urbanl * max_snow_layers))
    allocate(zc_improad_rb(num_urbanl * max_snow_layers))
    allocate(zc_perroad_rb(num_urbanl * max_snow_layers))

  end subroutine urbanxx_snowLayer_init

  !-----------------------------------------------------------------------
  subroutine PackActiveCount(num_urbanc, filter_urbanc, &
       buf_roof, buf_improad, buf_perroad)
    !
    ! !DESCRIPTION:
    ! Pack -snl(c) (number of active snow layers) for roof, impervious
    ! road, pervious road into flat per-surface buffers, indexed by a
    ! per-surface packed counter (same idx-increments-per-surface convention
    ! already used by e.g. SetImperviousRoadSoilWater in UrbanxxMod.F90).
    !
    use ColumnType    , only : col_pp
    use column_varcon , only : icol_roof, icol_road_imperv, icol_road_perv
    !
    implicit none
    integer(c_int), intent(in)  :: num_urbanc
    integer        , intent(in)  :: filter_urbanc(:)
    real(c_double) , intent(out) :: buf_roof(:), buf_improad(:), buf_perroad(:)

    integer :: fc, c
    integer :: idx_roof, idx_improad, idx_perroad

    idx_roof = 0
    idx_improad = 0
    idx_perroad = 0
    do fc = 1, num_urbanc
      c = filter_urbanc(fc)
      select case (col_pp%itype(c))
      case (icol_roof)
        idx_roof = idx_roof + 1
        buf_roof(idx_roof) = real(-col_pp%snl(c), c_double)
      case (icol_road_imperv)
        idx_improad = idx_improad + 1
        buf_improad(idx_improad) = real(-col_pp%snl(c), c_double)
      case (icol_road_perv)
        idx_perroad = idx_perroad + 1
        buf_perroad(idx_perroad) = real(-col_pp%snl(c), c_double)
      end select
    end do

  end subroutine PackActiveCount

  !-----------------------------------------------------------------------
  subroutine PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, &
       begc, elm_array, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    !
    ! !DESCRIPTION:
    ! Pack an ELM per-column, per-layer quantity (t_soisno, h2osoi_liq,
    ! h2osoi_ice, dz, or z) for the active snow-layer range of roof,
    ! impervious road, pervious road columns into flat per-surface buffers
    ! sized (num_urbanl * max_snow_layers), using the elm_index =
    ! -(cpp_index) mapping. Inactive slots (cpp_index >= -snl(c)) are
    ! zero-filled.
    !
    ! elm_array is declared with an explicit lower bound (begc, -nlevsno+1)
    ! rather than plain assumed-shape (:,:) — assumed-shape dummies default
    ! to 1-based indexing, which would silently break both the column index
    ! c (natively begc:endc, not 1:) and the negative snow-layer index
    ! elm_j. This mirrors the existing `tk( bounds%begc: , -nlevsno+1: )`
    ! declaration style already used in SoilTemperatureMod.F90.
    !
    use ColumnType    , only : col_pp
    use column_varcon , only : icol_roof, icol_road_imperv, icol_road_perv
    !
    implicit none
    integer(c_int), intent(in)  :: num_urbanc
    integer        , intent(in)  :: filter_urbanc(:)
    integer(c_int), intent(in)  :: num_urbanl
    integer        , intent(in)  :: begc
    real(r8)       , intent(in)  :: elm_array(begc:, -max_snow_layers+1:)
    logical        , intent(in)  :: isLayoutLeft
    real(c_double) , intent(out) :: buf_roof(:), buf_improad(:), buf_perroad(:)

    integer :: fc, c, j, n, elm_j, flat_idx
    integer :: idx_roof, idx_improad, idx_perroad
    real(r8) :: val

    idx_roof = 0
    idx_improad = 0
    idx_perroad = 0
    do fc = 1, num_urbanc
      c = filter_urbanc(fc)
      select case (col_pp%itype(c))
      case (icol_roof)
        idx_roof = idx_roof + 1
        call PackOneColumn(c, idx_roof, buf_roof)
      case (icol_road_imperv)
        idx_improad = idx_improad + 1
        call PackOneColumn(c, idx_improad, buf_improad)
      case (icol_road_perv)
        idx_perroad = idx_perroad + 1
        call PackOneColumn(c, idx_perroad, buf_perroad)
      end select
    end do

  contains

    subroutine PackOneColumn(col, idx, buf)
      integer, intent(in) :: col, idx
      real(c_double), intent(inout) :: buf(:)

      n = -col_pp%snl(col)
      do j = 0, max_snow_layers - 1
        if (j < n) then
          elm_j = -j
          val = elm_array(col, elm_j)
        else
          val = 0._r8
        end if
        if (isLayoutLeft) then
          flat_idx = j * num_urbanl + idx
        else
          flat_idx = (idx - 1) * max_snow_layers + j + 1
        end if
        buf(flat_idx) = real(val, c_double)
      end do
    end subroutine PackOneColumn

  end subroutine PackLayerQuantity

  !-----------------------------------------------------------------------
  subroutine PackLayerQuantityInt(num_urbanc, filter_urbanc, num_urbanl, &
       begc, elm_array, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    !
    ! !DESCRIPTION:
    ! Same as PackLayerQuantity but for an integer ELM per-layer array
    ! (e.g. col_ef%imelt), double-marshalled into the output buffers for
    ! the C API. Used by UrbanxxSnowCompactionMod.F90.
    !
    use ColumnType    , only : col_pp
    use column_varcon , only : icol_roof, icol_road_imperv, icol_road_perv
    !
    implicit none
    integer(c_int), intent(in)  :: num_urbanc
    integer        , intent(in)  :: filter_urbanc(:)
    integer(c_int), intent(in)  :: num_urbanl
    integer        , intent(in)  :: begc
    integer        , intent(in)  :: elm_array(begc:, -max_snow_layers+1:)
    logical        , intent(in)  :: isLayoutLeft
    real(c_double) , intent(out) :: buf_roof(:), buf_improad(:), buf_perroad(:)

    integer :: fc, c, j, n, elm_j, flat_idx
    integer :: idx_roof, idx_improad, idx_perroad
    integer :: val

    idx_roof = 0
    idx_improad = 0
    idx_perroad = 0
    do fc = 1, num_urbanc
      c = filter_urbanc(fc)
      select case (col_pp%itype(c))
      case (icol_roof)
        idx_roof = idx_roof + 1
        call PackOneColumnInt(c, idx_roof, buf_roof)
      case (icol_road_imperv)
        idx_improad = idx_improad + 1
        call PackOneColumnInt(c, idx_improad, buf_improad)
      case (icol_road_perv)
        idx_perroad = idx_perroad + 1
        call PackOneColumnInt(c, idx_perroad, buf_perroad)
      end select
    end do

  contains

    subroutine PackOneColumnInt(col, idx, buf)
      integer, intent(in) :: col, idx
      real(c_double), intent(inout) :: buf(:)

      n = -col_pp%snl(col)
      do j = 0, max_snow_layers - 1
        if (j < n) then
          elm_j = -j
          val = elm_array(col, elm_j)
        else
          val = 0
        end if
        if (isLayoutLeft) then
          flat_idx = j * num_urbanl + idx
        else
          flat_idx = (idx - 1) * max_snow_layers + j + 1
        end if
        buf(flat_idx) = real(val, c_double)
      end do
    end subroutine PackOneColumnInt

  end subroutine PackLayerQuantityInt

  !-----------------------------------------------------------------------
  subroutine SetAndReadBack1D(set_sub, get_sub, num_urbanl, buf, readback_buf)
    !
    ! !DESCRIPTION:
    ! Call a 1D UrbanSet* function then immediately the matching UrbanGet*
    ! function, storing the result in readback_buf for later comparison in
    ! urbanxx_snowLayer_check.
    !
    implicit none
    interface
      subroutine set_sub(urban, values, length, status)
        use iso_c_binding
        use urban_mod, only : UrbanType
        type(UrbanType), intent(in) :: urban
        type(c_ptr), value :: values
        integer(c_int), value :: length
        integer(c_int), intent(out) :: status
      end subroutine set_sub
      subroutine get_sub(urban, values, length, status)
        use iso_c_binding
        use urban_mod, only : UrbanType
        type(UrbanType), intent(in) :: urban
        type(c_ptr), value :: values
        integer(c_int), value :: length
        integer(c_int), intent(out) :: status
      end subroutine get_sub
    end interface
    integer(c_int), intent(in) :: num_urbanl
    real(c_double), target, intent(in)  :: buf(:)
    real(c_double), target, intent(inout) :: readback_buf(:)

    integer(c_int) :: status

    call set_sub(urbanxx, c_loc(buf), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call get_sub(urbanxx, c_loc(readback_buf), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine SetAndReadBack1D

  !-----------------------------------------------------------------------
  subroutine SetAndReadBack2D(set_sub, get_sub, num_urbanl, buf, readback_buf)
    !
    ! !DESCRIPTION:
    ! Same as SetAndReadBack1D but for the (num_urbanl, max_snow_layers) 2D
    ! quantities.
    !
    use elm_varpar , only : nlevsno
    !
    implicit none
    interface
      subroutine set_sub(urban, values, size, status)
        use iso_c_binding
        use urban_mod, only : UrbanType
        type(UrbanType), intent(in) :: urban
        type(c_ptr), value :: values
        integer(c_int), dimension(2), intent(in) :: size
        integer(c_int), intent(out) :: status
      end subroutine set_sub
      subroutine get_sub(urban, values, size, status)
        use iso_c_binding
        use urban_mod, only : UrbanType
        type(UrbanType), intent(in) :: urban
        type(c_ptr), value :: values
        integer(c_int), dimension(2), intent(in) :: size
        integer(c_int), intent(out) :: status
      end subroutine get_sub
    end interface
    integer(c_int), intent(in) :: num_urbanl
    real(c_double), target, intent(in)  :: buf(:)
    real(c_double), target, intent(inout) :: readback_buf(:)

    integer(c_int) :: status
    integer(c_int), dimension(2) :: size2D

    size2D(1) = num_urbanl
    size2D(2) = max_snow_layers

    call set_sub(urbanxx, c_loc(buf), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call get_sub(urbanxx, c_loc(readback_buf), size2D, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine SetAndReadBack2D

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowLayer(num_urbanl, num_urbanc, filter_urbanc, begc)
    !
    ! !DESCRIPTION:
    ! Pure plumbing round-trip: pack ELM's current active-count/thickness/Zc
    ! (all unaffected by the Phase 4 solve — see module-level comment), seed
    ! URBANxx, and immediately read it back into the round-trip buffers.
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
    real(c_double), allocatable, target :: buf_roof(:), buf_improad(:), buf_perroad(:)

    isLayoutLeft = UrbanKokkosIsLayoutLeft()

    allocate(buf_roof(num_urbanl))
    allocate(buf_improad(num_urbanl))
    allocate(buf_perroad(num_urbanl))
    call PackActiveCount(num_urbanc, filter_urbanc, buf_roof, buf_improad, buf_perroad)
    call SetAndReadBack1D(UrbanSetNumberOfActiveSnowLayersRoof, &
         UrbanGetNumberOfActiveSnowLayersRoof, num_urbanl, buf_roof, numactive_roof_rb)
    call SetAndReadBack1D(UrbanSetNumberOfActiveSnowLayersImperviousRoad, &
         UrbanGetNumberOfActiveSnowLayersImperviousRoad, num_urbanl, buf_improad, numactive_improad_rb)
    call SetAndReadBack1D(UrbanSetNumberOfActiveSnowLayersPerviousRoad, &
         UrbanGetNumberOfActiveSnowLayersPerviousRoad, num_urbanl, buf_perroad, numactive_perroad_rb)
    deallocate(buf_roof, buf_improad, buf_perroad)

    allocate(buf_roof(num_urbanl * max_snow_layers))
    allocate(buf_improad(num_urbanl * max_snow_layers))
    allocate(buf_perroad(num_urbanl * max_snow_layers))

    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_pp%dz, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call SetAndReadBack2D(UrbanSetSnowLayerThicknessRoof, UrbanGetSnowLayerThicknessRoof, &
         num_urbanl, buf_roof, dz_roof_rb)
    call SetAndReadBack2D(UrbanSetSnowLayerThicknessImperviousRoad, UrbanGetSnowLayerThicknessImperviousRoad, &
         num_urbanl, buf_improad, dz_improad_rb)
    call SetAndReadBack2D(UrbanSetSnowLayerThicknessPerviousRoad, UrbanGetSnowLayerThicknessPerviousRoad, &
         num_urbanl, buf_perroad, dz_perroad_rb)

    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_pp%z, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call SetAndReadBack2D(UrbanSetSnowLayerZcRoof, UrbanGetSnowLayerZcRoof, &
         num_urbanl, buf_roof, zc_roof_rb)
    call SetAndReadBack2D(UrbanSetSnowLayerZcImperviousRoad, UrbanGetSnowLayerZcImperviousRoad, &
         num_urbanl, buf_improad, zc_improad_rb)
    call SetAndReadBack2D(UrbanSetSnowLayerZcPerviousRoad, UrbanGetSnowLayerZcPerviousRoad, &
         num_urbanl, buf_perroad, zc_perroad_rb)

    deallocate(buf_roof, buf_improad, buf_perroad)

  end subroutine urbanxx_snowLayer

  !-----------------------------------------------------------------------
  subroutine CheckOne(label, fresh_buf, readback_buf, n)
    !
    ! !DESCRIPTION:
    ! Compare a freshly-packed ELM buffer against the value URBANxx handed
    ! back immediately after being set this same timestep. In Phase 1 there
    ! is no physics between the two, so this should match to bit precision;
    ! any nonzero difference means the Kokkos view layout/packing is wrong.
    !
    implicit none
    character(len=*), intent(in) :: label
    real(c_double)  , intent(in) :: fresh_buf(:)
    real(c_double)  , intent(in) :: readback_buf(:)
    integer          , intent(in) :: n

    integer :: i
    real(r8) :: max_err

    max_err = 0._r8
    do i = 1, n
      max_err = max(max_err, abs(fresh_buf(i) - readback_buf(i)))
    end do
    if (max_err > 0._r8) then
      write(iulog,'(A,A,A,ES24.16)') 'ERROR: UrbanxxSnowLayerMod round-trip mismatch for ', &
           trim(label), ', max_err = ', max_err
      call endrun(msg='UrbanxxSnowLayerMod: round-trip check failed for '//trim(label))
    end if

  end subroutine CheckOne

  !-----------------------------------------------------------------------
  subroutine urbanxx_snowLayer_check(num_urbanl, num_urbanc, filter_urbanc, begc)
    !
    ! !DESCRIPTION:
    ! Pure plumbing round-trip validation: repack ELM's (unchanged, since no
    ! physics ran in between) active-count/thickness/Zc and compare against
    ! the readback buffers populated by urbanxx_snowLayer this same
    ! timestep.
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
    real(c_double), allocatable, target :: buf_roof(:), buf_improad(:), buf_perroad(:)

    isLayoutLeft = UrbanKokkosIsLayoutLeft()

    allocate(buf_roof(num_urbanl))
    allocate(buf_improad(num_urbanl))
    allocate(buf_perroad(num_urbanl))
    call PackActiveCount(num_urbanc, filter_urbanc, buf_roof, buf_improad, buf_perroad)
    call CheckOne('NumberOfActiveSnowLayersRoof', buf_roof, numactive_roof_rb, num_urbanl)
    call CheckOne('NumberOfActiveSnowLayersImperviousRoad', buf_improad, numactive_improad_rb, num_urbanl)
    call CheckOne('NumberOfActiveSnowLayersPerviousRoad', buf_perroad, numactive_perroad_rb, num_urbanl)
    deallocate(buf_roof, buf_improad, buf_perroad)

    allocate(buf_roof(num_urbanl * max_snow_layers))
    allocate(buf_improad(num_urbanl * max_snow_layers))
    allocate(buf_perroad(num_urbanl * max_snow_layers))

    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_pp%dz, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call CheckOne('SnowLayerThicknessRoof', buf_roof, dz_roof_rb, num_urbanl * max_snow_layers)
    call CheckOne('SnowLayerThicknessImperviousRoad', buf_improad, dz_improad_rb, num_urbanl * max_snow_layers)
    call CheckOne('SnowLayerThicknessPerviousRoad', buf_perroad, dz_perroad_rb, num_urbanl * max_snow_layers)

    call PackLayerQuantity(num_urbanc, filter_urbanc, num_urbanl, begc, &
         col_pp%z, isLayoutLeft, buf_roof, buf_improad, buf_perroad)
    call CheckOne('SnowLayerZcRoof', buf_roof, zc_roof_rb, num_urbanl * max_snow_layers)
    call CheckOne('SnowLayerZcImperviousRoad', buf_improad, zc_improad_rb, num_urbanl * max_snow_layers)
    call CheckOne('SnowLayerZcPerviousRoad', buf_perroad, zc_perroad_rb, num_urbanl * max_snow_layers)

    deallocate(buf_roof, buf_improad, buf_perroad)

    if (masterproc) then
      write(iulog,*) 'UrbanxxSnowLayerMod: round-trip check passed for 3 quantities '// &
           '(NumberOfActiveSnowLayers, Thickness, Zc)'
    end if

  end subroutine urbanxx_snowLayer_check

end module UrbanxxSnowLayerMod
