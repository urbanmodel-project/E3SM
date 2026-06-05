module UrbanxxSnowMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Snow-related URBANxx driver routines.
  !-----------------------------------------------------------------------
  use iso_c_binding
  use urban_mod
  use abortutils          , only : endrun
  use elm_varctl          , only : iulog
  use UrbanxxInstanceMod  , only : urbanxx

  implicit none
  private

  public :: urbanxx_updateSnowFraction

contains

  subroutine urbanxx_updateSnowFraction(num_urbanl, filter_urbanl)
    integer, intent(in) :: num_urbanl
    integer, intent(in) :: filter_urbanl(:)
    !
    integer(c_int) :: status
    !
    if (num_urbanl == 0) return

    call UrbanComputeUpdateSnowFraction(urbanxx, status)
    if (status /= URBAN_SUCCESS) then
       write(iulog,*) 'urbanxx_updateSnowFraction: UrbanComputeUpdateSnowFraction failed, status=', status
       call endrun('urbanxx_updateSnowFraction: URBANxx error')
    end if

  end subroutine urbanxx_updateSnowFraction

end module UrbanxxSnowMod
