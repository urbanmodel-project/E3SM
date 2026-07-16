module UrbanxxEnergyBalanceMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Energy-balance-closure (errsoi) diagnostic for URBANxx. Mirrors ELM's
  ! own self-consistency check (Phasechange_beta's xmf(c) accumulator +
  ! SoilFluxesMod's errsoi_patch residual + BalanceCheckMod's tolerance
  ! check), but computed entirely from URBANxx's OWN state: its own
  ! accumulated phase-change energy (Xmf), its own pre/post-solve layer
  ! temperatures, and its own EflxSoilGrnd. This validates that URBANxx's
  ! own coupled temperature-diffusion + phase-change solve conserves
  ! energy — it does NOT compare against ELM's xmf/errsoi_col. Per
  ! plans/CLAUDE.md, ELM and URBANxx are independently-implemented,
  ! intentionally non-BFB solvers that are allowed to diverge; this check
  ! answers a different question than the ELM-vs-URBANxx value comparisons
  ! in the other Urbanxx*_check routines ("is URBANxx's own solve
  ! internally energy-conserving?"), not "do ELM and URBANxx agree?".
  !
  ! v1 scope (see plans/ELM_URBANxx_Known_Gaps_And_Future_Work.md, "No
  ! xmf-equivalent energy-balance-closure diagnostic"):
  !   - All 5 surfaces: roof, sunlit/shaded wall, impervious/pervious road.
  !     Walls carry no water/ice state (WallDataType has no water-content
  !     views at all), so Xmf is identically 0 for walls, hardcoded in the
  !     residual formula below rather than plumbed through a dead C++ view.
  !   - Roof/impervious-road/pervious-road columns with active discrete
  !     snow layers (NumberOfActiveSnowLayers > 0) are skipped: the
  !     snow-layer sensible-heat-storage term needs a per-layer heat
  !     capacity that is recomputed transiently inside
  !     ApplyMultiLayerSnowPhaseChange from wliq/wice and never persisted,
  !     so a faithful recompute would need further C++-side plumbing.
  !     Skipping (not approximating) avoids a check that silently reports a
  !     spurious residual on every snowy timestep. Walls have no snow
  !     layers at all, so no such skip is needed for sunwall/shadewall.
  !   - Roof/walls use ComputeHeatDiffusion's building-surface boundary
  !     condition (useTopLayerAdjustment=false), so unlike roads their
  !     "fact" needs no dz_eff/Zc/Zi geometry — it is simply
  !     dtime/cv_times_dz for every layer including the top one (confirmed
  !     against UrbanHeatDiffusion.cpp's Solve1DHeatDiffusion/
  !     ApplyThinSnowPhaseChange, both called with useTopLayerAdjustment=
  !     false for roof/sunwall/shadewall). Roof still needs the same
  !     transient thin-snow cv_times_dz(l,0) correction addback as
  !     impervious road (see h2osno_roof below) since that correction is
  !     independent of the dz_eff geometry.
  !   - EflxBuildingHeat (roof/sunwall/shadewall only, per ELM's own
  !     formula which scopes that term to icol_sunwall/icol_shadewall/
  !     icol_roof) is now persisted C++-side (previously a discarded local
  !     in UrbanHeatDiffusion.cpp).
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use shr_kind_mod               , only : r8 => shr_kind_r8
  use spmdMod                    , only : masterproc, iam
  use elm_varctl                 , only : iulog
  use UrbanxxInstanceMod         , only : urbanxx
  use UrbanxxSoilFluxesMod       , only : eflx_soil_grnd_improad, eflx_soil_grnd_perroad, &
                                          eflx_soil_grnd_roof, eflx_soil_grnd_sunwall, &
                                          eflx_soil_grnd_shadwall
  use UrbanxxSoilTemperatureMod  , only : layertemp_improad, layertemp_perroad, &
                                          layertemp_roof, layertemp_sunwall, layertemp_shadwall
  use UrbanxxDiagnosticControlMod, only : urbanxx_softfail_permits

  implicit none

  private

  ! Persistent buffers (allocated once in init)
  real(c_double) , allocatable, target :: xmf_improad(:), xmf_perroad(:), xmf_roof(:)
  real(c_double) , allocatable, target :: told_improad(:), told_perroad(:)
  real(c_double) , allocatable, target :: told_roof(:), told_sunwall(:), told_shadwall(:)
  real(c_double) , allocatable, target :: cvtimesdz_improad(:), cvtimesdz_perroad(:)
  real(c_double) , allocatable, target :: cvtimesdz_roof(:), cvtimesdz_sunwall(:), cvtimesdz_shadwall(:)
  real(c_double) , allocatable, target :: dz_improad(:), dz_perroad(:)
  real(c_double) , allocatable, target :: zc_improad(:), zc_perroad(:)
  real(c_double) , allocatable, target :: zi_improad(:), zi_perroad(:)
  real(c_double) , allocatable, target :: numactivesnow_improad(:), numactivesnow_perroad(:)
  real(c_double) , allocatable, target :: numactivesnow_roof(:)
  ! Roof/impervious-road-only: needed to reproduce the transient thin-snow
  ! heat-capacity correction to cv_times_dz(l,0) that ComputeHeatDiffusion
  ! applies then reverts within the same timestep — see the getter's doc
  ! comment in UrbanHeatDiffusionGetters.cpp for why pervious road needs no
  ! equivalent (its own correction is never reverted within the timestep).
  real(c_double) , allocatable, target :: h2osno_improad(:), h2osno_roof(:)
  ! Building-interior heat flux term (roof/sunwall/shadewall only).
  real(c_double) , allocatable, target :: eflxbldgheat_roof(:)
  real(c_double) , allocatable, target :: eflxbldgheat_sunwall(:)
  real(c_double) , allocatable, target :: eflxbldgheat_shadwall(:)

  public :: urbanxx_energyBalance_init
  public :: urbanxx_energyBalance
  public :: urbanxx_energyBalance_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_energyBalance_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for the energy-balance-closure check.
    ! Called once during initialization.
    !
    use elm_varpar, only : nlevgrnd, nlevurb
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl

    allocate(xmf_improad(num_urbanl))
    allocate(xmf_perroad(num_urbanl))
    allocate(xmf_roof(num_urbanl))
    allocate(told_improad(num_urbanl * nlevgrnd))
    allocate(told_perroad(num_urbanl * nlevgrnd))
    allocate(told_roof(num_urbanl * nlevurb))
    allocate(told_sunwall(num_urbanl * nlevurb))
    allocate(told_shadwall(num_urbanl * nlevurb))
    allocate(cvtimesdz_improad(num_urbanl * nlevgrnd))
    allocate(cvtimesdz_perroad(num_urbanl * nlevgrnd))
    allocate(cvtimesdz_roof(num_urbanl * nlevurb))
    allocate(cvtimesdz_sunwall(num_urbanl * nlevurb))
    allocate(cvtimesdz_shadwall(num_urbanl * nlevurb))
    allocate(dz_improad(num_urbanl * nlevgrnd))
    allocate(dz_perroad(num_urbanl * nlevgrnd))
    allocate(zc_improad(num_urbanl * nlevgrnd))
    allocate(zc_perroad(num_urbanl * nlevgrnd))
    allocate(zi_improad(num_urbanl * (nlevgrnd + 1)))
    allocate(zi_perroad(num_urbanl * (nlevgrnd + 1)))
    allocate(numactivesnow_improad(num_urbanl))
    allocate(numactivesnow_perroad(num_urbanl))
    allocate(numactivesnow_roof(num_urbanl))
    allocate(h2osno_improad(num_urbanl))
    allocate(h2osno_roof(num_urbanl))
    allocate(eflxbldgheat_roof(num_urbanl))
    allocate(eflxbldgheat_sunwall(num_urbanl))
    allocate(eflxbldgheat_shadwall(num_urbanl))

  end subroutine urbanxx_energyBalance_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_energyBalance(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Fetch URBANxx's own Xmf, pre-solve temperature, and layer-geometry
    ! state needed by urbanxx_energyBalance_check. Must be called after
    ! urbanxx_soilFluxes (so EflxSoilGrnd, read directly from
    ! UrbanxxSoilFluxesMod, is fresh) and after urbanxx_soilTemperature (so
    ! Xmf/TOld and layertemp_improad/perroad, the post-solve temperature,
    ! reflect this timestep's solve).
    !
    use elm_varpar, only : nlevgrnd, nlevurb
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    !
    integer(c_int)               :: status
    integer(c_int), dimension(2) :: size2D_soil, size2D_zi, size2D_urban

    size2D_soil(1) = num_urbanl
    size2D_soil(2) = nlevgrnd
    size2D_zi(1)   = num_urbanl
    size2D_zi(2)   = nlevgrnd + 1
    size2D_urban(1) = num_urbanl
    size2D_urban(2) = nlevurb

    call UrbanGetXmfImperviousRoad(urbanxx, c_loc(xmf_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetXmfPerviousRoad(urbanxx, c_loc(xmf_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetXmfRoof(urbanxx, c_loc(xmf_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetTemperatureOldImperviousRoad(urbanxx, c_loc(told_improad), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetTemperatureOldPerviousRoad(urbanxx, c_loc(told_perroad), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetTemperatureOldRoof(urbanxx, c_loc(told_roof), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetTemperatureOldSunlitWall(urbanxx, c_loc(told_sunwall), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetTemperatureOldShadedWall(urbanxx, c_loc(told_shadwall), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetCvTimesDzImperviousRoad(urbanxx, c_loc(cvtimesdz_improad), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCvTimesDzPerviousRoad(urbanxx, c_loc(cvtimesdz_perroad), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCvTimesDzRoof(urbanxx, c_loc(cvtimesdz_roof), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCvTimesDzSunlitWall(urbanxx, c_loc(cvtimesdz_sunwall), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCvTimesDzShadedWall(urbanxx, c_loc(cvtimesdz_shadwall), size2D_urban, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetDzImperviousRoad(urbanxx, c_loc(dz_improad), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetDzPerviousRoad(urbanxx, c_loc(dz_perroad), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetZcImperviousRoad(urbanxx, c_loc(zc_improad), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetZcPerviousRoad(urbanxx, c_loc(zc_perroad), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetZiImperviousRoad(urbanxx, c_loc(zi_improad), size2D_zi, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetZiPerviousRoad(urbanxx, c_loc(zi_perroad), size2D_zi, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetNumberOfActiveSnowLayersImperviousRoad(urbanxx, c_loc(numactivesnow_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetNumberOfActiveSnowLayersPerviousRoad(urbanxx, c_loc(numactivesnow_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetNumberOfActiveSnowLayersRoof(urbanxx, c_loc(numactivesnow_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetH2OSnoImperviousRoad(urbanxx, c_loc(h2osno_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetH2OSnoRoof(urbanxx, c_loc(h2osno_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetEflxBuildingHeatRoof(urbanxx, c_loc(eflxbldgheat_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetEflxBuildingHeatSunlitWall(urbanxx, c_loc(eflxbldgheat_sunwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetEflxBuildingHeatShadedWall(urbanxx, c_loc(eflxbldgheat_shadwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine urbanxx_energyBalance

  !-----------------------------------------------------------------------
  integer function FlatIndex(l, j, n1, n2, isLeft) result(idx)
    !
    ! !DESCRIPTION:
    ! 1-based flat index into a (n1 x n2) array returned by a 2D URBANxx
    ! getter, matching the packing convention used throughout the
    ! Urbanxx*Mod files (e.g. UrbanxxSoilTemperatureMod's
    ! layertemp_improad unflattening): l varies fastest under Kokkos
    ! LayoutLeft, j varies fastest under LayoutRight.
    !
    implicit none
    integer, intent(in) :: l, j, n1, n2   ! 1-based: l in [1,n1], j in [1,n2]
    logical, intent(in) :: isLeft

    if (isLeft) then
      idx = l + (j - 1) * n1
    else
      idx = j + (l - 1) * n2
    end if

  end function FlatIndex

  !-----------------------------------------------------------------------
  subroutine urbanxx_energyBalance_check(num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Compute the errsoi residual for all 5 urban surfaces (roof,
    ! sunlit/shaded wall, impervious/pervious road) — roof/impervious/
    ! pervious-road columns with active discrete snow layers are skipped
    ! (v1 scope, see module header) — and check it against a fixed physical
    ! tolerance (not an ELM value) using the same urbanxx_softfail_permits
    ! ceiling used by every other Urbanxx*_check routine.
    !
    use ColumnType             , only : col_pp
    use column_varcon          , only : icol_road_imperv, icol_road_perv, &
                                        icol_roof, icol_sunwall, icol_shadewall
    use elm_varpar              , only : nlevgrnd, nlevurb
    use timeInfoMod            , only : dtime_mod
    use elm_time_manager        , only : get_nstep, get_curr_date
    use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
    use elm_varcon              , only : cpice
    !
    implicit none
    !
    integer(c_int), intent(in) :: num_urbanl
    integer       , intent(in) :: filter_urbanl(:)
    integer(c_int), intent(in) :: num_urbanc
    integer       , intent(in) :: filter_urbanc(:)
    !
    real(r8), parameter :: capr     = 0.34_r8    ! matches UrbanHeatDiffusion.cpp's Turing factor
    real(r8), parameter :: tol_warn  = 1.0e-5_r8  ! W/m^2 — starting guess, matches ELM's own errsoi warn threshold
    real(r8), parameter :: tol_abort = 1.0e-4_r8  ! W/m^2 — starting guess, matches ELM's own errsoi hard-abort threshold
    !
    logical  :: isLayoutLeft
    integer  :: fc, c, j, l, idx_improad, idx_perroad, idx_roof, idx_sunwall, idx_shadewall
    real(r8) :: errsoi, layer_sum, fact, dz2, dz_eff
    real(r8) :: max_error
    real(r8) :: eflx_val, xmf_val, tnew, told, cv, dz0, dz1_geom, zi0

    ! Cold-start guard: TOld/Xmf are not meaningful before the first real
    ! ComputeHeatDiffusion solve has run with fully-seeded state — same
    ! nstep()<=2 rationale as urbanxx_netShortwave_check.
    if (get_nstep() <= 2) return

    isLayoutLeft = UrbanKokkosIsLayoutLeft()
    max_error = 0._r8
    idx_improad = 0
    idx_perroad = 0
    idx_roof = 0
    idx_sunwall = 0
    idx_shadewall = 0

    do fc = 1, num_urbanc
      c = filter_urbanc(fc)

      select case (col_pp%itype(c))
      case (icol_road_imperv)
        idx_improad = idx_improad + 1
        l = idx_improad
        if (nint(numactivesnow_improad(l)) /= 0) cycle  ! v1: no active-snow coverage yet

        eflx_val = eflx_soil_grnd_improad(l)
        xmf_val  = xmf_improad(l)
        layer_sum = 0._r8
        do j = 1, nlevgrnd
          tnew = layertemp_improad(FlatIndex(l, j, num_urbanl, nlevgrnd, isLayoutLeft))
          told = told_improad(FlatIndex(l, j, num_urbanl, nlevgrnd, isLayoutLeft))
          cv   = cvtimesdz_improad(FlatIndex(l, j, num_urbanl, nlevgrnd, isLayoutLeft))
          if (j == 1) then
            ! Reproduce the transient thin-snow heat-capacity correction
            ! ComputeHeatDiffusion applies to layer 0 then reverts within
            ! the same timestep (Step 3c in UrbanHeatDiffusion.cpp) — the
            ! CvTimesDzImperviousRoad getter reads the POST-reversion
            ! value, so it must be added back here to match what the
            ! solve/phase-change actually used. N==0 is already guaranteed
            ! at this point (cycled above otherwise).
            if (h2osno_improad(l) > 0._r8) cv = cv + cpice * h2osno_improad(l)
            dz0      = dz_improad(FlatIndex(l, 1, num_urbanl, nlevgrnd, isLayoutLeft))
            zi0      = zi_improad(FlatIndex(l, 1, num_urbanl, nlevgrnd + 1, isLayoutLeft))
            dz1_geom = zc_improad(FlatIndex(l, 1, num_urbanl, nlevgrnd, isLayoutLeft)) - zi0
            dz2      = zc_improad(FlatIndex(l, 2, num_urbanl, nlevgrnd, isLayoutLeft)) - zi0
            dz_eff   = 0.5_r8 * (dz1_geom + capr * dz2)
            fact     = dtime_mod / cv * dz0 / dz_eff
          else
            fact = dtime_mod / cv
          end if
          layer_sum = layer_sum + (tnew - told) / fact
        end do
        errsoi = eflx_val - xmf_val - layer_sum
        call CheckErrsoi('ImperviousRoad', c, errsoi, tol_warn, tol_abort, max_error)

      case (icol_road_perv)
        idx_perroad = idx_perroad + 1
        l = idx_perroad
        if (nint(numactivesnow_perroad(l)) /= 0) cycle  ! v1: no active-snow coverage yet

        eflx_val = eflx_soil_grnd_perroad(l)
        xmf_val  = xmf_perroad(l)
        layer_sum = 0._r8
        do j = 1, nlevgrnd
          tnew = layertemp_perroad(FlatIndex(l, j, num_urbanl, nlevgrnd, isLayoutLeft))
          told = told_perroad(FlatIndex(l, j, num_urbanl, nlevgrnd, isLayoutLeft))
          cv   = cvtimesdz_perroad(FlatIndex(l, j, num_urbanl, nlevgrnd, isLayoutLeft))
          if (j == 1) then
            dz0      = dz_perroad(FlatIndex(l, 1, num_urbanl, nlevgrnd, isLayoutLeft))
            zi0      = zi_perroad(FlatIndex(l, 1, num_urbanl, nlevgrnd + 1, isLayoutLeft))
            dz1_geom = zc_perroad(FlatIndex(l, 1, num_urbanl, nlevgrnd, isLayoutLeft)) - zi0
            dz2      = zc_perroad(FlatIndex(l, 2, num_urbanl, nlevgrnd, isLayoutLeft)) - zi0
            dz_eff   = 0.5_r8 * (dz1_geom + capr * dz2)
            fact     = dtime_mod / cv * dz0 / dz_eff
          else
            fact = dtime_mod / cv
          end if
          layer_sum = layer_sum + (tnew - told) / fact
        end do
        errsoi = eflx_val - xmf_val - layer_sum
        call CheckErrsoi('PerviousRoad', c, errsoi, tol_warn, tol_abort, max_error)

      case (icol_roof)
        idx_roof = idx_roof + 1
        l = idx_roof
        if (nint(numactivesnow_roof(l)) /= 0) cycle  ! v1: no active-snow coverage yet

        eflx_val = eflx_soil_grnd_roof(l)
        xmf_val  = xmf_roof(l)
        layer_sum = 0._r8
        do j = 1, nlevurb
          tnew = layertemp_roof(FlatIndex(l, j, num_urbanl, nlevurb, isLayoutLeft))
          told = told_roof(FlatIndex(l, j, num_urbanl, nlevurb, isLayoutLeft))
          cv   = cvtimesdz_roof(FlatIndex(l, j, num_urbanl, nlevurb, isLayoutLeft))
          ! Roof is a building surface (useTopLayerAdjustment=false in
          ! UrbanHeatDiffusion.cpp), so unlike roads there is no dz_eff
          ! adjustment for the top layer — fact = dtime/cv uniformly. The
          ! transient thin-snow cv_times_dz(l,0) correction still needs its
          ! addback though (independent of the dz_eff geometry); N==0 is
          ! already guaranteed at this point (cycled above otherwise).
          if (j == 1 .and. h2osno_roof(l) > 0._r8) cv = cv + cpice * h2osno_roof(l)
          fact = dtime_mod / cv
          layer_sum = layer_sum + (tnew - told) / fact
        end do
        errsoi = eflx_val - xmf_val + eflxbldgheat_roof(l) - layer_sum
        call CheckErrsoi('Roof', c, errsoi, tol_warn, tol_abort, max_error)

      case (icol_sunwall)
        idx_sunwall = idx_sunwall + 1
        l = idx_sunwall

        eflx_val = eflx_soil_grnd_sunwall(l)
        layer_sum = 0._r8
        do j = 1, nlevurb
          tnew = layertemp_sunwall(FlatIndex(l, j, num_urbanl, nlevurb, isLayoutLeft))
          told = told_sunwall(FlatIndex(l, j, num_urbanl, nlevurb, isLayoutLeft))
          cv   = cvtimesdz_sunwall(FlatIndex(l, j, num_urbanl, nlevurb, isLayoutLeft))
          fact = dtime_mod / cv
          layer_sum = layer_sum + (tnew - told) / fact
        end do
        ! Xmf = 0: walls carry no water/ice state, so phase change there is
        ! physically impossible (per plans/CLAUDE.md-documented decision —
        ! no dead Kokkos view added for it).
        errsoi = eflx_val + eflxbldgheat_sunwall(l) - layer_sum
        call CheckErrsoi('SunlitWall', c, errsoi, tol_warn, tol_abort, max_error)

      case (icol_shadewall)
        idx_shadewall = idx_shadewall + 1
        l = idx_shadewall

        eflx_val = eflx_soil_grnd_shadwall(l)
        layer_sum = 0._r8
        do j = 1, nlevurb
          tnew = layertemp_shadwall(FlatIndex(l, j, num_urbanl, nlevurb, isLayoutLeft))
          told = told_shadwall(FlatIndex(l, j, num_urbanl, nlevurb, isLayoutLeft))
          cv   = cvtimesdz_shadwall(FlatIndex(l, j, num_urbanl, nlevurb, isLayoutLeft))
          fact = dtime_mod / cv
          layer_sum = layer_sum + (tnew - told) / fact
        end do
        errsoi = eflx_val + eflxbldgheat_shadwall(l) - layer_sum  ! Xmf = 0, see SunlitWall
        call CheckErrsoi('ShadedWall', c, errsoi, tol_warn, tol_abort, max_error)

      end select
    end do

    write(iulog,*) 'Max error in energy balance (errsoi, W/m^2): ', max_error

  contains

    !-----------------------------------------------------------------------
    subroutine CheckErrsoi(surfname, c, errsoi, tol_warn, tol_abort, max_error)
      !
      ! !DESCRIPTION:
      ! Shared warn/abort logic for a single computed errsoi residual,
      ! following the same soft-fail-wrapped pattern as every other
      ! Urbanxx*_check routine.
      !
      implicit none
      character(len=*), intent(in)    :: surfname
      integer          , intent(in)    :: c
      real(r8)         , intent(in)    :: errsoi, tol_warn, tol_abort
      real(r8)         , intent(inout) :: max_error
      !
      integer :: yr, mo, dy, tod

      max_error = max(max_error, abs(errsoi))

      if (abs(errsoi) > tol_warn) then
        call get_curr_date(yr, mo, dy, tod)
        write(iulog,'(A)') 'WARNING: URBANxx energy-balance-closure (errsoi) residual exceeds tolerance'
        write(iulog,'(A,A)')     '  surface       = ', trim(surfname)
        write(iulog,'(A,I0)')    '  column        = ', c
        write(iulog,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '  step=', get_nstep(), ', date=', yr, '-', mo, '-', dy, ', tod=', tod
        write(iulog,'(A,ES24.16)') '  errsoi (W/m^2) = ', errsoi

        if (abs(errsoi) > tol_abort) then
          if (urbanxx_softfail_permits(abs(errsoi), tol_abort)) then
            write(iulog,'(A)') '  URBANxx check failure permitted to continue (urbanxx_softfail/urbanxx_always_continue); logging and continuing'
          else
            call exit(0)
          end if
        end if
      end if

    end subroutine CheckErrsoi

  end subroutine urbanxx_energyBalance_check

end module UrbanxxEnergyBalanceMod
