module UrbanxxSoilTemperatureMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Set building temperature and compute heat diffusion
  ! for urban areas in the Urban++ model.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use LandunitDataType     , only : lun_es
  use UrbanxxInstanceMod   , only : urbanxx

  implicit none

  private

  ! Persistent input buffers (allocated once in init)
  real(c_double) , allocatable, target :: buildingTemp(:)

  ! Persistent output buffers for layer temperatures (allocated once in init)
  ! Roof, sunlit wall, shaded wall: (num_urbanl * nlevurb) — 1D flat arrays
  ! Impervious road, pervious road: (num_urbanl * nlevgrnd) — 1D flat arrays
  real(c_double) , allocatable, target, public :: layertemp_roof(:)
  real(c_double) , allocatable, target, public :: layertemp_improad(:)
  real(c_double) , allocatable, target, public :: layertemp_perroad(:)
  real(c_double) , allocatable, target, public :: layertemp_sunwall(:)
  real(c_double) , allocatable, target, public :: layertemp_shadwall(:)

  public :: urbanxx_soilTemperature_init
  public :: urbanxx_soilTemperature

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilTemperature_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for heat diffusion computation.
    ! Called once during initialization.
    !
    use elm_varpar, only : nlevgrnd, nlevurb
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl

    ! Input buffers
    allocate(buildingTemp(num_urbanl))

    ! Output buffers — 2D flattened
    allocate(layertemp_roof(num_urbanl * nlevurb))
    allocate(layertemp_sunwall(num_urbanl * nlevurb))
    allocate(layertemp_shadwall(num_urbanl * nlevurb))
    allocate(layertemp_improad(num_urbanl * nlevgrnd))
    allocate(layertemp_perroad(num_urbanl * nlevgrnd))

  end subroutine urbanxx_soilTemperature_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilTemperature(num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Set building temperature and compute heat diffusion for urban areas
    !
    use TemperatureType, only : temperature_type
    use ColumnType     , only : col_pp
    use column_varcon  , only : icol_road_perv
    use elm_varpar     , only : nlevgrnd, nlevurb
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int)         , intent(in) :: num_urbanl        ! number of urban landunits
    integer(c_int)         , intent(in) :: num_urbanc        ! number of urban columns
    integer                , intent(in) :: filter_urbanl(:)  ! urban layer filter
    integer                , intent(in) :: filter_urbanc(:)  ! urban column filter
    type(temperature_type) , intent(in) :: temperature_vars
    !
    ! !LOCAL VARIABLES:
    integer(c_int)                       :: status
    integer                              :: fc, c, j, fl, l
    integer(c_int), dimension(2)         :: size2D_urban, size2D_soil

    ! Set building temperature before heat diffusion
    associate(                          &
        t_building => lun_es%t_building & ! Input: [real(r8) (:)   ]  internal building temperature (K)
        )

      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         buildingTemp(fl) = t_building(l)
      end do

      call UrbanSetBuildingTemperature(urbanxx, c_loc(buildingTemp), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanComputeHeatDiffusion(urbanxx, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! Extract layer temperatures from UrbanXX
      ! Roof, sunlit wall, shaded wall: (num_urbanl, nlevurb)
      size2D_urban(1) = num_urbanl
      size2D_urban(2) = nlevurb

      call UrbanGetLayerTempRoof(urbanxx, c_loc(layertemp_roof), size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetLayerTempSunlitWall(urbanxx, c_loc(layertemp_sunwall), size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetLayerTempShadedWall(urbanxx, c_loc(layertemp_shadwall), size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! Impervious road, pervious road: (num_urbanl, nlevgrnd)
      size2D_soil(1) = num_urbanl
      size2D_soil(2) = nlevgrnd

      call UrbanGetLayerTempImperviousRoad(urbanxx, c_loc(layertemp_improad), size2D_soil, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetLayerTempPerviousRoad(urbanxx, c_loc(layertemp_perroad), size2D_soil, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    end associate

  end subroutine urbanxx_soilTemperature

end module UrbanxxSoilTemperatureMod
