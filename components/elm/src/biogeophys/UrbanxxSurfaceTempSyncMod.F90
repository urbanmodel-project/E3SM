module UrbanxxSurfaceTempSyncMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Discrete multi-layer snowpack: Phase 4, surface-temperature sync fix for
  ! Known Gap #1. See plans/ELM_URBANxx_MultiLayerSnow_Port_Scoping_Plan.md.
  !
  ! Root cause this closes: ELM creates a new discrete snow layer in
  ! CanopyHydrologyMod.F90 (snl 0 -> -1, elm_driver.F90 step 1), well before
  ! urbanxx_surfaceFluxes runs (step 2). ELM's own t_grnd(c) — the value its
  ! surface-flux formulas use — is recomputed in CanopyTemperatureMod.F90
  ! (called between CanopyHydrology and urbanxx_surfaceFluxes) as a
  ! frac_sno_eff-weighted blend of the new snow layer's temperature and the
  ! top soil layer's temperature. URBANxx's EffectiveSurfTemp, in contrast,
  ! is only ever updated inside ComputeHeatDiffusion (urbanxx_soilTemperature,
  ! step 3) — so at the moment urbanxx_surfaceFluxes runs, it still reflects
  ! the PREVIOUS timestep's post-solve value, with no knowledge of a snow
  ! layer ELM created a few driver-steps earlier in the SAME timestep. This
  ! is exactly the mismatch that trips the urbanxx_surfaceFluxes_check
  ! tolerance abort (see the scoping doc's "Why this doc exists now").
  !
  ! Fix: a new call site (elm_driver.F90, between CanopyTemperature and
  ! UrbanFluxes/urbanxx_surfaceFluxes) syncs EffectiveSurfTemp to ELM's
  ! freshly-recomputed t_grnd(c) for roof, impervious road, and pervious
  ! road (the 3 surfaces that accumulate snow) every timestep — mirroring
  ! the same "temporary, over-seed" philosophy Phase 1 already established
  ! for the snow-layer state itself (see UrbanxxSnowLayerMod.F90), just
  ! applied earlier in the timestep and to this one additional scalar.
  ! TGrnd0 does not need the same treatment: it is only consumed inside
  ! UrbanHeatDiffusion.cpp (after ComputeHeatDiffusion re-derives it from
  ! current state), never by urbanxx_surfaceFluxes.
  !
  ! Second call site (Phase 6, added once seasonal-cycle validation found
  ! the urbanxx_netLongwave check failing at the first timestep after real
  ! snow-layer restructuring): elm_driver.F90 calls this routine AGAIN right
  ! after HydrologyNoDrainage, which is where ELM performs its OWN last
  ! t_grnd(c) update for the timestep (HydrologyNoDrainageMod.F90, after
  ! CombineSnowLayers/DivideSnowLayers may have merged/split snow layers).
  ! Without this second call, EffectiveSurfTemp would carry the
  ! PRE-restructuring value set by ComputeHeatDiffusion (step 3) forward into
  ! the NEXT timestep's urbanxx_netLongwave — which runs before the first
  ! call site (after CanopyTemperature) fires again — reading a stale
  ! topmost-layer temperature ELM has already moved on from. Both call sites
  ! push the same thing (col_es%t_grnd for roof/impervious/pervious road);
  ! this routine doesn't need to know which point in the timestep it's being
  ! called from.
  !
  ! Walls (sunlit/shaded) are deliberately excluded: they never accumulate
  ! snow (see the scoping doc's resolved "which surfaces accumulate
  ! multi-layer snow" question), so their EffectiveSurfTemp is unaffected by
  ! this issue.

  use iso_c_binding
  use urban_mod
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use spmdMod            , only : iam
  use UrbanxxInstanceMod , only : urbanxx

  implicit none

  private

  public :: urbanxx_syncSurfaceTemp

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_syncSurfaceTemp(num_urbanl, num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Pack ELM's current t_grnd(c) for roof/impervious road/pervious road
    ! columns and push it into URBANxx's EffectiveSurfTemp, before
    ! urbanxx_surfaceFluxes reads it this same timestep.
    !
    use ColumnDataType , only : col_es
    use ColumnType      , only : col_pp
    use column_varcon   , only : icol_roof, icol_road_imperv, icol_road_perv
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer        , intent(in) :: filter_urbanc(:)

    real(c_double), allocatable, target :: buf_roof(:), buf_improad(:), buf_perroad(:)
    integer(c_int) :: status
    integer :: fc, c
    integer :: idx_roof, idx_improad, idx_perroad

    allocate(buf_roof(num_urbanl))
    allocate(buf_improad(num_urbanl))
    allocate(buf_perroad(num_urbanl))

    idx_roof = 0
    idx_improad = 0
    idx_perroad = 0
    do fc = 1, num_urbanc
      c = filter_urbanc(fc)
      select case (col_pp%itype(c))
      case (icol_roof)
        idx_roof = idx_roof + 1
        buf_roof(idx_roof) = real(col_es%t_grnd(c), c_double)
      case (icol_road_imperv)
        idx_improad = idx_improad + 1
        buf_improad(idx_improad) = real(col_es%t_grnd(c), c_double)
      case (icol_road_perv)
        idx_perroad = idx_perroad + 1
        buf_perroad(idx_perroad) = real(col_es%t_grnd(c), c_double)
      end select
    end do

    call UrbanSetEffectiveSurfaceTempRoof(urbanxx, c_loc(buf_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetEffectiveSurfaceTempImperviousRoad(urbanxx, c_loc(buf_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanSetEffectiveSurfaceTempPerviousRoad(urbanxx, c_loc(buf_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    deallocate(buf_roof, buf_improad, buf_perroad)

  end subroutine urbanxx_syncSurfaceTemp

end module UrbanxxSurfaceTempSyncMod
