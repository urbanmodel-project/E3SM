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

  end subroutine urbanxx_surfaceFluxes_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_surfaceFluxes(num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, &
       surfalb_vars, urbanparams_vars, frictionvel_vars)
    !
    implicit none
    !
    integer(c_int)     , intent(in) :: num_urbanl
    integer            , intent(in) :: filter_urbanl(:)         ! urban landunit filter
    integer(c_int)     , intent(in) :: num_urbanc
    integer            , intent(in) :: filter_urbanc(:)         ! urban column filter
    type(surfalb_type) , intent(in) :: surfalb_vars
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(frictionvel_type) , intent(in)    :: frictionvel_vars
    !
    integer(c_int)                       :: status

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
