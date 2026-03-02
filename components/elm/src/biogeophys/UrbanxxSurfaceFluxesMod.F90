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
  use ColumnDataType       , only : col_ws
  use UrbanxxInstanceMod   , only : urbanxx
  use UrbanxxMod           , only : SetHeightParameters
  use VegetationDataType   , only : veg_ef, veg_wf
  use VegetationType       , only : veg_pp
  use LandunitDataType     , only : lun_es, lun_ws

  implicit none

  private

  ! Persistent input buffers (allocated once in init)
  real(c_double) , allocatable, target :: fwet_road(:)
  real(c_double) , allocatable, target :: fwet_roof(:)

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

    ! Input buffers
    allocate(fwet_road(num_urbanl))
    allocate(fwet_roof(num_urbanl))

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

  end subroutine urbanxx_surfaceFluxes_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_surfaceFluxes(num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, &
       num_urbanp, filter_urbanp, surfalb_vars, urbanparams_vars, frictionvel_vars)
    !
    use column_varcon, only : icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, icol_shadewall
    use ColumnType, only : col_pp
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
    real(r8)                             :: max_error
    real(r8)                             :: err_eflx_sh_grnd, err_qflx_evap_soi, err_cgrnds, err_cgrndl
    integer :: idx_roof, idx_road_improv, idx_road_perv, idx_sunwall, idx_shadwall, idx_landunit

    call SetHeightParameters(urbanxx, num_urbanl, filter_urbanl, &
       urbanparams_vars, frictionvel_vars)
    call SetFwetValues(urbanxx, num_urbanl, num_urbanc, filter_urbanc)

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
      case (icol_road_imperv)
         idx_road_improv = idx_road_improv + 1
         err_eflx_sh_grnd  = abs(veg_ef%eflx_sh_grnd(p)   - eflx_sh_grnd_improad(idx_road_improv))
         err_qflx_evap_soi = abs(veg_wf%qflx_evap_soi(p)  - qflx_evap_soi_improad(idx_road_improv))
         err_cgrnds        = abs(veg_ef%cgrnds(p)          - cgrnds_improad(idx_road_improv))
         err_cgrndl        = abs(veg_ef%cgrndl(p)          - cgrndl_improad(idx_road_improv))
      case (icol_road_perv)
         idx_road_perv = idx_road_perv + 1
         err_eflx_sh_grnd  = abs(veg_ef%eflx_sh_grnd(p)   - eflx_sh_grnd_perroad(idx_road_perv))
         err_qflx_evap_soi = abs(veg_wf%qflx_evap_soi(p)  - qflx_evap_soi_perroad(idx_road_perv))
         err_cgrnds        = abs(veg_ef%cgrnds(p)          - cgrnds_perroad(idx_road_perv))
         err_cgrndl        = abs(veg_ef%cgrndl(p)          - cgrndl_perroad(idx_road_perv))
      case (icol_sunwall)
         idx_sunwall = idx_sunwall + 1
         err_eflx_sh_grnd = abs(veg_ef%eflx_sh_grnd(p)    - eflx_sh_grnd_sunwall(idx_sunwall))
         err_cgrnds       = abs(veg_ef%cgrnds(p)           - cgrnds_sunwall(idx_sunwall))
         err_cgrndl       = abs(veg_ef%cgrndl(p)           - cgrndl_sunwall(idx_sunwall))
      case (icol_shadewall)
         idx_shadwall = idx_shadwall + 1
         err_eflx_sh_grnd = abs(veg_ef%eflx_sh_grnd(p)    - eflx_sh_grnd_shadwall(idx_shadwall))
         err_cgrnds       = abs(veg_ef%cgrnds(p)           - cgrnds_shadwall(idx_shadwall))
         err_cgrndl       = abs(veg_ef%cgrndl(p)           - cgrndl_shadwall(idx_shadwall))
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
   write(iulog,*) 'Max error in SH ground      : ', err_eflx_sh_grnd
   write(iulog,*) 'Max error in Evap soil      : ', err_qflx_evap_soi
   write(iulog,*) 'Max error in d(SH)/dT       : ', err_cgrnds
   write(iulog,*) 'Max error in d(Evap)/dT     : ', err_cgrndl

   max_error = 0._r8
   do fl = 1, num_urbanl
      l = filter_urbanl(fl)
      max_error = max(max_error, abs(lun_es%taf(l) - taf(fl)))
      max_error = max(max_error, abs(lun_ws%qaf(l) - qaf(fl)))
   end do
   write(iulog,*) 'Max error in surface Taf, Qaf       : ', max_error
   if (max_error > 1.0e-10) then
      write(iulog,*) 'Error exceeds tolerance! Check Urban++ surface Taf/Qaf computation.'
      call exit(0)
   end if
  end subroutine urbanxx_surfaceFluxes

  !-----------------------------------------------------------------------
  subroutine SetFwetValues(urban, num_urbanl, num_urbanc, filter_urbanc)
    !
    ! Set fraction wet values for impervious road based on snow depth and ponding
    !
    use WaterStateType, only : waterstate_type
    use elm_varcon, only : pondmx_urban
    use ColumnType, only : col_pp
    use column_varcon, only : icol_roof, icol_road_imperv
    !
    implicit none
    !
    type(UrbanType)      , intent(in) :: urban
    integer(c_int)       , intent(in) :: num_urbanl
    integer(c_int)       , intent(in) :: num_urbanc
    integer              , intent(in) :: filter_urbanc(:) ! urban column filter
    !
    integer(c_int)                       :: status
    integer                              :: fc, l, c, idx_road, idx_roof
    real(r8)                             :: fwet

    associate(                             &
         snow_depth => col_ws%snow_depth , & ! Input: [real(r8) (:)] snow depth (m)
         h2osoi_liq => col_ws%h2osoi_liq , & ! Input: [real(r8) (:,:)] liquid water (kg/m2)
         h2osoi_ice => col_ws%h2osoi_ice   & ! Input: [real(r8) (:,:)] ice water (kg/m2
         )

      ! Loop through urban landunits
      idx_road = 0
      idx_roof = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)

        if (col_pp%itype(c) == icol_roof) then
           idx_roof = idx_roof + 1
           if (snow_depth(c) > 0._r8) then
              fwet = min(snow_depth(c)/0.05_r8, 1._r8)
           else
              fwet = (max(0._r8, h2osoi_liq(c,1)+h2osoi_ice(c,1))/pondmx_urban)**0.666666666666_r8
              fwet = min(fwet,1._r8)
           end if
           fwet_roof(idx_roof) = fwet
        end if
        if (col_pp%itype(c) == icol_road_imperv) then
           idx_road = idx_road + 1
           ! Calculate fraction wet based on snow depth or ponding
           ! From UrbanFluxesMod.F90:L577-582
           if (snow_depth(c) > 0._r8) then
              fwet = min(snow_depth(c)/0.05_r8, 1._r8)
           else
              fwet = (max(0._r8, h2osoi_liq(c,1)+h2osoi_ice(c,1))/pondmx_urban)**0.666666666666_r8
              fwet = min(fwet, 1._r8)
          end if
          fwet_road(idx_road) = fwet
        endif
      end do

      ! Set the values in urbanxx
      call UrbanSetFractionWetImperviousRoad(urban, c_loc(fwet_road), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetFractionWetRoof(urban, c_loc(fwet_roof), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    end associate

  end subroutine SetFwetValues

end module UrbanxxSurfaceFluxesMod
