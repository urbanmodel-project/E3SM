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

  public :: urbanxx_netShortwave

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_netShortwave(num_urbanl, filter_urbanl, surfalb_vars, &
         urbanparams_vars, frictionvel_vars)
    !
    implicit none
    !
    integer(c_int)         , intent(in)  :: num_urbanl
    integer                , intent(in)  :: filter_urbanl(:)         ! urban landunit filter
    type(surfalb_type)     , intent(in)  :: surfalb_vars
    type(urbanparams_type) , intent(in)  :: urbanparams_vars
    type(frictionvel_type) , intent(in)  :: frictionvel_vars
    !
    integer(c_int)                       :: status
    integer                              :: fl, l
    real(c_double) , allocatable, target :: atmCoszen(:)

    associate(                  &
         coli =>    lun_pp%coli & ! Input:  [integer (:)    ]  beginning column index for landunit
         )

      allocate(atmCoszen(num_urbanl))
      ! Fill arrays with values from ELM data structures
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         atmCoszen(fl)   = surfalb_vars%coszen_col(coli(l))  ! Assumes coszen for each column are the same
      end do

      call UrbanSetAtmCoszen(urbanxx, c_loc(atmCoszen), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanComputeNetShortwave(urbanxx, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      deallocate(atmCoszen)

    end associate
  end subroutine urbanxx_netShortwave

end module UrbanxxNetShortwaveMod
