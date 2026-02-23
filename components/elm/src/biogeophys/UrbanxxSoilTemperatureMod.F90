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

  public :: urbanxx_soilTemperature

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilTemperature(num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Placeholder for soil temperature calculations in urban areas
    !
    use TemperatureType, only : temperature_type
    use ColumnType     , only : col_pp
    use column_varcon  , only : icol_road_perv
    use elm_varpar     , only : nlevgrnd
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
    real(c_double) , allocatable, target :: buildingTemp(:)

    ! Set building temperature before heat diffusion
    associate(                          &
        t_building => lun_es%t_building & ! Input: [real(r8) (:)   ]  internal building temperature (K)
        )

      allocate(buildingTemp(num_urbanl))
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         buildingTemp(fl) = t_building(l)
      end do

      call UrbanSetBuildingTemperature(urbanxx, c_loc(buildingTemp), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      deallocate(buildingTemp)

      call UrbanComputeHeatDiffusion(urbanxx, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    end associate

  end subroutine urbanxx_soilTemperature

end module UrbanxxSoilTemperatureMod
