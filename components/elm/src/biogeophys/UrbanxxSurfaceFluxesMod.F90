module UrbanxxSurfaceFluxesMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute surface fluxes for the Urban++ model.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use urban_kokkos_interface
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use UrbanParamsType      , only : urbanparams_type
  use SurfaceAlbedoType    , only : surfalb_type
  use FrictionVelocityType , only : frictionvel_type
  use UrbanxxInstanceMod   , only : urbanxx
  use UrbanxxMod           , only : SetHeightParameters
  use VegetationDataType   , only : veg_ef, veg_wf
  use VegetationType       , only : veg_pp
  use LandunitDataType     , only : lun_es, lun_ws

  implicit none

  private

  ! Persistent output buffers (allocated once in init)
  real(c_double) , allocatable, target, public :: eflx_sh_grnd_roof(:)
  real(c_double) , allocatable, target, public :: eflx_sh_grnd_improad(:)
  real(c_double) , allocatable, target, public :: eflx_sh_grnd_perroad(:)
  real(c_double) , allocatable, target, public :: eflx_sh_grnd_sunwall(:)
  real(c_double) , allocatable, target, public :: eflx_sh_grnd_shadwall(:)
  real(c_double) , allocatable, target, public :: qflx_evap_soi_roof(:)
  real(c_double) , allocatable, target, public :: qflx_evap_soi_improad(:)
  real(c_double) , allocatable, target, public :: qflx_evap_soi_perroad(:)
  real(c_double) , allocatable, target, public :: cgrnds_roof(:)
  real(c_double) , allocatable, target, public :: cgrnds_improad(:)
  real(c_double) , allocatable, target, public :: cgrnds_perroad(:)
  real(c_double) , allocatable, target, public :: cgrnds_sunwall(:)
  real(c_double) , allocatable, target, public :: cgrnds_shadwall(:)
  real(c_double) , allocatable, target, public :: cgrndl_roof(:)
  real(c_double) , allocatable, target, public :: cgrndl_improad(:)
  real(c_double) , allocatable, target, public :: cgrndl_perroad(:)
  real(c_double) , allocatable, target, public :: cgrndl_sunwall(:)
  real(c_double) , allocatable, target, public :: cgrndl_shadwall(:)
  real(c_double) , allocatable, target, public :: taf(:)
  real(c_double) , allocatable, target, public :: qaf(:)

  ! Persistent input buffers: TopH2OSoiLiq/Ice synced from ELM before each surface flux call
  real(c_double) , allocatable, target :: in_top_liq_roof(:)
  real(c_double) , allocatable, target :: in_top_ice_roof(:)
  real(c_double) , allocatable, target :: in_top_liq_improad(:)
  real(c_double) , allocatable, target :: in_top_ice_improad(:)

  public :: urbanxx_surfaceFluxes_init
  public :: urbanxx_surfaceFluxes

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_surfaceFluxes_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for surface fluxes computation.
    ! Called once during initialization.
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl

    ! Output buffers
    allocate(eflx_sh_grnd_roof(num_urbanl))
    allocate(eflx_sh_grnd_improad(num_urbanl))
    allocate(eflx_sh_grnd_perroad(num_urbanl))
    allocate(eflx_sh_grnd_sunwall(num_urbanl))
    allocate(eflx_sh_grnd_shadwall(num_urbanl))
    allocate(qflx_evap_soi_roof(num_urbanl))
    allocate(qflx_evap_soi_improad(num_urbanl))
    allocate(qflx_evap_soi_perroad(num_urbanl))
    allocate(cgrnds_roof(num_urbanl))
    allocate(cgrnds_improad(num_urbanl))
    allocate(cgrnds_perroad(num_urbanl))
    allocate(cgrnds_sunwall(num_urbanl))
    allocate(cgrnds_shadwall(num_urbanl))
    allocate(cgrndl_roof(num_urbanl))
    allocate(cgrndl_improad(num_urbanl))
    allocate(cgrndl_perroad(num_urbanl))
    allocate(cgrndl_sunwall(num_urbanl))
    allocate(cgrndl_shadwall(num_urbanl))
    allocate(taf(num_urbanl))
    allocate(qaf(num_urbanl))

    ! Input buffers for TopH2OSoiLiq/Ice sync
    allocate(in_top_liq_roof(num_urbanl))
    allocate(in_top_ice_roof(num_urbanl))
    allocate(in_top_liq_improad(num_urbanl))
    allocate(in_top_ice_improad(num_urbanl))

  end subroutine urbanxx_surfaceFluxes_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_surfaceFluxes(num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, &
       num_urbanp, filter_urbanp, surfalb_vars, urbanparams_vars, frictionvel_vars)
    !
    use column_varcon, only : icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, icol_shadewall
    use ColumnType, only : col_pp
    use ColumnDataType, only : col_ws
    !
    implicit none
    !
    integer(c_int)     , intent(in) :: num_urbanl
    integer            , intent(in) :: filter_urbanl(:)         ! urban landunit filter
    integer(c_int)     , intent(in) :: num_urbanc
    integer            , intent(in) :: filter_urbanc(:)         ! urban column filter
    integer(c_int)     , intent(in) :: num_urbanp
    integer            , intent(in) :: filter_urbanp(:)         ! urban point filter
    type(surfalb_type) , intent(in) :: surfalb_vars
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(frictionvel_type) , intent(in)    :: frictionvel_vars
    !
    integer(c_int)                       :: status
    integer                              :: fp, c, p, fl, l
    real(r8)                             :: max_error, max_rel_error
    real(r8)                             :: err_eflx_sh_grnd, err_qflx_evap_soi, err_cgrnds, err_cgrndl
    real(r8)                             :: rel_eflx_sh_grnd, rel_qflx_evap_soi, rel_cgrnds, rel_cgrndl
    integer :: idx_roof, idx_road_improv, idx_road_perv, idx_sunwall, idx_shadwall, idx_landunit
    integer :: fc
    integer :: idx_top_roof, idx_top_improad

    ! --- Sync TopH2OSoiLiq/Ice from ELM state before surface flux computation ---
    ! URBANxx's stored TopH2OSoiLiq/Ice may be stale (from urbanxx_dewCondensation)
    ! while ELM's h2osoi_liq(c,1) has been further updated by WaterTable (dew addition).
    ! Syncing here ensures fwet is computed from the correct water state, which
    ! in turn gives the correct QflxEvapSoil that matches ELM's qflx_evap_soi.
    idx_top_roof    = 0
    idx_top_improad = 0
    do fc = 1, num_urbanc
      c = filter_urbanc(fc)
      select case (col_pp%itype(c))
      case (icol_roof)
        idx_top_roof = idx_top_roof + 1
        in_top_liq_roof(idx_top_roof)    = real(col_ws%h2osoi_liq(c,1), c_double)
        in_top_ice_roof(idx_top_roof)    = real(col_ws%h2osoi_ice(c,1), c_double)
      case (icol_road_imperv)
        idx_top_improad = idx_top_improad + 1
        in_top_liq_improad(idx_top_improad) = real(col_ws%h2osoi_liq(c,1), c_double)
        in_top_ice_improad(idx_top_improad) = real(col_ws%h2osoi_ice(c,1), c_double)
      end select
    end do

    !call UrbanSetTopH2OSoiLiqRoof(urbanxx, c_loc(in_top_liq_roof), &
    !     num_urbanl, status)
    !if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    !call UrbanSetTopH2OSoiIceRoof(urbanxx, c_loc(in_top_ice_roof), &
    !     num_urbanl, status)
    !if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    !call UrbanSetTopH2OSoiLiqImperviousRoad(urbanxx, c_loc(in_top_liq_improad), &
    !     num_urbanl, status)
    !if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    !call UrbanSetTopH2OSoiIceImperviousRoad(urbanxx, c_loc(in_top_ice_improad), &
    !     num_urbanl, status)
    !if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanComputeSurfaceFluxes(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Extract sensible heat flux from UrbanXX
    call UrbanGetSensibleHeatFluxRoof(urbanxx, c_loc(eflx_sh_grnd_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSensibleHeatFluxImperviousRoad(urbanxx, c_loc(eflx_sh_grnd_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSensibleHeatFluxPerviousRoad(urbanxx, c_loc(eflx_sh_grnd_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSensibleHeatFluxSunlitWall(urbanxx, c_loc(eflx_sh_grnd_sunwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetSensibleHeatFluxShadedWall(urbanxx, c_loc(eflx_sh_grnd_shadwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Extract soil evaporation flux from UrbanXX
    call UrbanGetEvapFluxRoof(urbanxx, c_loc(qflx_evap_soi_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetEvapFluxImperviousRoad(urbanxx, c_loc(qflx_evap_soi_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetEvapFluxPerviousRoad(urbanxx, c_loc(qflx_evap_soi_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Extract Cgrnds (d(sensible heat flux)/dT) from UrbanXX
    call UrbanGetCgrndsRoof(urbanxx, c_loc(cgrnds_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCgrndsImperviousRoad(urbanxx, c_loc(cgrnds_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCgrndsPerviousRoad(urbanxx, c_loc(cgrnds_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCgrndsSunlitWall(urbanxx, c_loc(cgrnds_sunwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCgrndsShadedWall(urbanxx, c_loc(cgrnds_shadwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Extract Cgrndl (d(latent heat flux)/dT) from UrbanXX
    call UrbanGetCgrndlRoof(urbanxx, c_loc(cgrndl_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCgrndlImperviousRoad(urbanxx, c_loc(cgrndl_improad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCgrndlPerviousRoad(urbanxx, c_loc(cgrndl_perroad), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCgrndlSunlitWall(urbanxx, c_loc(cgrndl_sunwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCgrndlShadedWall(urbanxx, c_loc(cgrndl_shadwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Extract canyon air temperature and humidity from UrbanXX
    call UrbanGetCanyonAirTemperature(urbanxx, c_loc(taf), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    call UrbanGetCanyonAirHumidity(urbanxx, c_loc(qaf), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

   idx_roof = 0
   idx_road_improv = 0
   idx_road_perv = 0
   idx_sunwall = 0
   idx_shadwall = 0

   err_eflx_sh_grnd  = 0._r8
   err_qflx_evap_soi = 0._r8
   err_cgrnds        = 0._r8
   err_cgrndl        = 0._r8
   rel_eflx_sh_grnd  = 0._r8
   rel_qflx_evap_soi = 0._r8
   rel_cgrnds        = 0._r8
   rel_cgrndl        = 0._r8

   do fp = 1, num_urbanp

      p = filter_urbanp(fp)
      c = veg_pp%column(p)

      select case (col_pp%itype(c))
      case (icol_roof)
         idx_roof = idx_roof + 1
         err_eflx_sh_grnd  = abs(veg_ef%eflx_sh_grnd(p)   - eflx_sh_grnd_roof(idx_roof))
         err_qflx_evap_soi = abs(veg_wf%qflx_evap_soi(p)  - qflx_evap_soi_roof(idx_roof))
         err_cgrnds        = abs(veg_ef%cgrnds(p)          - cgrnds_roof(idx_roof))
         err_cgrndl        = abs(veg_ef%cgrndl(p)          - cgrndl_roof(idx_roof))
         rel_eflx_sh_grnd  = err_eflx_sh_grnd  / max(abs(veg_ef%eflx_sh_grnd(p)),  1.0e-20_r8)
         rel_qflx_evap_soi = err_qflx_evap_soi / max(abs(veg_wf%qflx_evap_soi(p)), 1.0e-20_r8)
         rel_cgrnds        = err_cgrnds        / max(abs(veg_ef%cgrnds(p)),          1.0e-20_r8)
         rel_cgrndl        = err_cgrndl        / max(abs(veg_ef%cgrndl(p)),          1.0e-20_r8)
      case (icol_road_imperv)
         idx_road_improv = idx_road_improv + 1
         err_eflx_sh_grnd  = abs(veg_ef%eflx_sh_grnd(p)   - eflx_sh_grnd_improad(idx_road_improv))
         err_qflx_evap_soi = abs(veg_wf%qflx_evap_soi(p)  - qflx_evap_soi_improad(idx_road_improv))
         err_cgrnds        = abs(veg_ef%cgrnds(p)          - cgrnds_improad(idx_road_improv))
         err_cgrndl        = abs(veg_ef%cgrndl(p)          - cgrndl_improad(idx_road_improv))
         rel_eflx_sh_grnd  = err_eflx_sh_grnd  / max(abs(veg_ef%eflx_sh_grnd(p)),  1.0e-20_r8)
         rel_qflx_evap_soi = err_qflx_evap_soi / max(abs(veg_wf%qflx_evap_soi(p)), 1.0e-20_r8)
         rel_cgrnds        = err_cgrnds        / max(abs(veg_ef%cgrnds(p)),          1.0e-20_r8)
         rel_cgrndl        = err_cgrndl        / max(abs(veg_ef%cgrndl(p)),          1.0e-20_r8)
      case (icol_road_perv)
         idx_road_perv = idx_road_perv + 1
         err_eflx_sh_grnd  = abs(veg_ef%eflx_sh_grnd(p)   - eflx_sh_grnd_perroad(idx_road_perv))
         err_qflx_evap_soi = abs(veg_wf%qflx_evap_soi(p)  - qflx_evap_soi_perroad(idx_road_perv))
         err_cgrnds        = abs(veg_ef%cgrnds(p)          - cgrnds_perroad(idx_road_perv))
         err_cgrndl        = abs(veg_ef%cgrndl(p)          - cgrndl_perroad(idx_road_perv))
         rel_eflx_sh_grnd  = err_eflx_sh_grnd  / max(abs(veg_ef%eflx_sh_grnd(p)),  1.0e-20_r8)
         rel_qflx_evap_soi = err_qflx_evap_soi / max(abs(veg_wf%qflx_evap_soi(p)), 1.0e-20_r8)
         rel_cgrnds        = err_cgrnds        / max(abs(veg_ef%cgrnds(p)),          1.0e-20_r8)
         rel_cgrndl        = err_cgrndl        / max(abs(veg_ef%cgrndl(p)),          1.0e-20_r8)
      case (icol_sunwall)
         idx_sunwall = idx_sunwall + 1
         err_eflx_sh_grnd = abs(veg_ef%eflx_sh_grnd(p)    - eflx_sh_grnd_sunwall(idx_sunwall))
         err_cgrnds       = abs(veg_ef%cgrnds(p)           - cgrnds_sunwall(idx_sunwall))
         err_cgrndl       = abs(veg_ef%cgrndl(p)           - cgrndl_sunwall(idx_sunwall))
         rel_eflx_sh_grnd = err_eflx_sh_grnd / max(abs(veg_ef%eflx_sh_grnd(p)), 1.0e-20_r8)
         rel_cgrnds       = err_cgrnds       / max(abs(veg_ef%cgrnds(p)),        1.0e-20_r8)
         rel_cgrndl       = err_cgrndl       / max(abs(veg_ef%cgrndl(p)),        1.0e-20_r8)
      case (icol_shadewall)
         idx_shadwall = idx_shadwall + 1
         err_eflx_sh_grnd = abs(veg_ef%eflx_sh_grnd(p)    - eflx_sh_grnd_shadwall(idx_shadwall))
         err_cgrnds       = abs(veg_ef%cgrnds(p)           - cgrnds_shadwall(idx_shadwall))
         err_cgrndl       = abs(veg_ef%cgrndl(p)           - cgrndl_shadwall(idx_shadwall))
         rel_eflx_sh_grnd = err_eflx_sh_grnd / max(abs(veg_ef%eflx_sh_grnd(p)), 1.0e-20_r8)
         rel_cgrnds       = err_cgrnds       / max(abs(veg_ef%cgrnds(p)),        1.0e-20_r8)
         rel_cgrndl       = err_cgrndl       / max(abs(veg_ef%cgrndl(p)),        1.0e-20_r8)
      end select

      if (err_eflx_sh_grnd > 1.0e-10_r8 .or. err_qflx_evap_soi > 1.0e-10_r8 .or. &
          err_cgrnds       > 1.0e-10_r8 .or. err_cgrndl        > 1.0e-10_r8) then
         write(iulog,*) 'ERROR: Surface flux error exceeds tolerance at p=', p, ' c=', c
         write(iulog,*) '  err_eflx_sh_grnd  = ', err_eflx_sh_grnd
         write(iulog,*) '  err_qflx_evap_soi = ', err_qflx_evap_soi
         write(iulog,*) '  err_cgrnds        = ', err_cgrnds
         write(iulog,*) '  err_cgrndl        = ', err_cgrndl
         call exit(0)
      end if

   end do
   write(iulog,*) 'Max error in SH ground              : ', err_eflx_sh_grnd,  ' (rel: ', rel_eflx_sh_grnd,  ')'
   write(iulog,*) 'Max error in Evap soil              : ', err_qflx_evap_soi, ' (rel: ', rel_qflx_evap_soi, ')'
   write(iulog,*) 'Max error in d(SH)/dT               : ', err_cgrnds,        ' (rel: ', rel_cgrnds,        ')'
   write(iulog,*) 'Max error in d(Evap)/dT             : ', err_cgrndl,        ' (rel: ', rel_cgrndl,        ')'

   max_error = 0._r8
   max_rel_error = 0._r8
   do fl = 1, num_urbanl
      l = filter_urbanl(fl)
      max_error     = max(max_error,     abs(lun_es%taf(l) - taf(fl)))
      max_error     = max(max_error,     abs(lun_ws%qaf(l) - qaf(fl)))
      max_rel_error = max(max_rel_error, abs(lun_es%taf(l) - taf(fl)) / max(abs(lun_es%taf(l)), 1.0e-20_r8))
      max_rel_error = max(max_rel_error, abs(lun_ws%qaf(l) - qaf(fl)) / max(abs(lun_ws%qaf(l)), 1.0e-20_r8))
   end do
   write(iulog,*) 'Max error in surface Taf, Qaf       : ', max_error,    ' (rel: ', max_rel_error, ')'
   if (max_error > 1.0e-10) then
      write(iulog,*) 'Error exceeds tolerance! Check Urban++ surface Taf/Qaf computation.'
      call exit(0)
   end if
  end subroutine urbanxx_surfaceFluxes

end module UrbanxxSurfaceFluxesMod
