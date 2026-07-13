module UrbanxxNetShortwaveMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute net shortwave radiation for the Urban++ model.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use UrbanParamsType      , only : urbanparams_type
  use SurfaceAlbedoType    , only : surfalb_type
  use LandunitType         , only : lun_pp
  use FrictionVelocityType , only : frictionvel_type
  use UrbanxxInstanceMod   , only : urbanxx

  implicit none

  private

  ! Persistent input buffers (allocated once in init)
  real(c_double) , allocatable, target :: atmCoszen(:)

  ! Persistent output buffers (allocated once in init)
  real(c_double) , allocatable, target, public :: swnet_roof(:)
  real(c_double) , allocatable, target, public :: swnet_improad(:)
  real(c_double) , allocatable, target, public :: swnet_perroad(:)
  real(c_double) , allocatable, target, public :: swnet_sunwall(:)
  real(c_double) , allocatable, target, public :: swnet_shadwall(:)

  public :: urbanxx_netShortwave_init
  public :: urbanxx_netShortwave
  public :: urbanxx_netShortwave_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_netShortwave_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for net shortwave computation.
    ! Called once during initialization.
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl

    ! Input buffers
    allocate(atmCoszen(num_urbanl))

    ! Output buffers
    allocate(swnet_roof(num_urbanl))
    allocate(swnet_improad(num_urbanl))
    allocate(swnet_perroad(num_urbanl))
    allocate(swnet_sunwall(num_urbanl))
    allocate(swnet_shadwall(num_urbanl))

  end subroutine urbanxx_netShortwave_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_netShortwave(num_urbanl, filter_urbanl)
    !
    implicit none
    !
    integer(c_int)         , intent(in)  :: num_urbanl
    integer                , intent(in)  :: filter_urbanl(:)         ! urban landunit filter
    !
    integer(c_int)                       :: status


    call UrbanComputeNetShortwave(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Extract net shortwave radiation from UrbanXX
    call UrbanGetNetShortwaveRoof(urbanxx, c_loc(swnet_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetNetShortwaveImperviousRoad(urbanxx, c_loc(swnet_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetNetShortwavePerviousRoad(urbanxx, c_loc(swnet_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetNetShortwaveSunlitWall(urbanxx, c_loc(swnet_sunwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetNetShortwaveShadedWall(urbanxx, c_loc(swnet_shadwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine urbanxx_netShortwave

  !-----------------------------------------------------------------------
  subroutine urbanxx_netShortwave_check(num_urbanp, filter_urbanp, solarabs_vars)
    !
    ! !DESCRIPTION:
    ! Diagnostic-only comparison of URBANxx net shortwave radiation against
    ! ELM's per-surface absorbed ground shortwave (solarabs_vars%sabg_patch(p),
    ! set by UrbanRadiation and consumed by SoilTemperatureMod as
    ! sabg_lyr_col(c,lyr_top) for urban columns). Note: col_ef%eflx_sabg_lyr is
    ! NOT usable here — it is populated only under the PETSc external thermal
    ! model and only for non-urban columns, so it stays at spval for all urban
    ! columns. Does not modify ELM or URBANxx state.
    !
    use ColumnType     , only : col_pp
    use VegetationType , only : veg_pp
    use column_varcon  , only : icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, icol_shadewall
    use SolarAbsorbedType , only : solarabs_type
    use elm_time_manager , only : get_nstep, get_curr_date
    use elm_varcon        , only : spval
    !
    implicit none
    !
    integer(c_int)       , intent(in) :: num_urbanp
    integer              , intent(in) :: filter_urbanp(:)
    type(solarabs_type)  , intent(in) :: solarabs_vars
    !
    integer  :: fp, p, c, idx_landunit
    integer  :: idx_roof, idx_road_imperv, idx_road_perv, idx_sunwall, idx_shadwall
    real(r8) :: max_error, max_rel_error
    real(r8) :: elm_value, urbanxx_value, abs_diff, rel_diff, signed_diff
    integer  :: year, month, day, tod

    ! This routine runs inside elm_driver's "!$OMP PARALLEL DO" clump loop, so
    ! any module-level "have we run once yet" flag would be racy across
    ! clumps/threads. doalb (which gates every caller of this routine) is a
    ! single value shared by all clumps for a given elm_drv call, so instead
    ! key the skip on get_nstep() itself: before UrbanAlbedo's first-ever
    ! doalb cycle, ELM's sabs_*_dir/dif absorption coefficients are still
    ! their InitCold value (0, or pre-Init spval for some landunits), so
    ! sabg(p) for that cycle is not a meaningful quantity to compare against
    ! on either side (URBANxx's own AbsorbedShortRad is equally not yet
    ! meaningful before ComputeNetShortwave has ever run). get_nstep() is a
    ! pure query with no shared mutable state, so this is race-free.
    if (get_nstep() <= 2) return

    associate( &
         sabg => solarabs_vars%sabg_patch & ! Input: [real(r8) (:)] solar radiation absorbed by ground (W/m**2)
         )

      max_error       = 0._r8
      max_rel_error   = 0._r8
      idx_roof        = 0
      idx_road_imperv = 0
      idx_road_perv   = 0
      idx_sunwall     = 0
      idx_shadwall    = 0

      do fp = 1, num_urbanp
        p = filter_urbanp(fp)
        c = veg_pp%column(p)
        idx_landunit = col_pp%landunit(c)

        select case (col_pp%itype(c))
        case (icol_roof)
          idx_roof = idx_roof + 1
          elm_value     = sabg(p)
          urbanxx_value = swnet_roof(idx_roof)
        case (icol_sunwall)
          idx_sunwall = idx_sunwall + 1
          elm_value     = sabg(p)
          urbanxx_value = swnet_sunwall(idx_sunwall)
        case (icol_shadewall)
          idx_shadwall = idx_shadwall + 1
          elm_value     = sabg(p)
          urbanxx_value = swnet_shadwall(idx_shadwall)
        case (icol_road_imperv)
          idx_road_imperv = idx_road_imperv + 1
          elm_value     = sabg(p)
          urbanxx_value = swnet_improad(idx_road_imperv)
        case (icol_road_perv)
          idx_road_perv = idx_road_perv + 1
          elm_value     = sabg(p)
          urbanxx_value = swnet_perroad(idx_road_perv)
        case default
          cycle
        end select

        ! On cold start (and until the first doalb cycle has run), ELM's
        ! sabs_*_dir/dif absorption coefficients are still spval, so sabg(p)
        ! is not yet a meaningful value to compare against. Skip rather than
        ! flag a spurious mismatch.
        if (elm_value == spval) cycle

        abs_diff     = abs(elm_value - urbanxx_value)
        rel_diff     = abs_diff / max(abs(elm_value), 1.0e-20_r8)
        signed_diff  = elm_value - urbanxx_value
        max_error     = max(max_error, abs_diff)
        max_rel_error = max(max_rel_error, rel_diff)

        if (abs_diff > 1.0e-9_r8) then
          call get_curr_date(year, month, day, tod)
          write(iulog,'(A)') 'ERROR: URBANxx net-shortwave comparison exceeds tolerance'
          write(iulog,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '  step=', get_nstep(), ', date=', year, '-', month, '-', day, ', tod=', tod
          write(iulog,'(A,I0,A,I0,A,I0,A,I0)') '  patch=', p, ', column=', c, ', landunit=', idx_landunit, ', column_type=', col_pp%itype(c)
          write(iulog,'(A,ES24.16)') '  ELM sabg(p)            =', elm_value
          write(iulog,'(A,ES24.16)') '  URBANxx net shortwave  =', urbanxx_value
          write(iulog,'(A,ES24.16)') '  signed difference      =', signed_diff
          write(iulog,'(A,ES24.16)') '  absolute difference    =', abs_diff
          write(iulog,'(A,ES24.16)') '  relative difference    =', rel_diff
          call exit(0)
        end if
      end do

      write(iulog,*) 'Max error in net shortwave          : ', max_error, ' (rel: ', max_rel_error, ')'

    end associate

  end subroutine urbanxx_netShortwave_check

end module UrbanxxNetShortwaveMod
