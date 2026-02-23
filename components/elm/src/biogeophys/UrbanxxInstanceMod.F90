module UrbanxxInstanceMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Holds the shared UrbanType instance and radiation constants used by
  ! all Urbanxx physics modules.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod, only : UrbanType

  implicit none

  private

  ! Constants for radiation bands and types
  integer(c_int), parameter, public :: numBands = 2  ! VIS, NIR
  integer(c_int), parameter, public :: numTypes = 2  ! Direct, Diffuse

  ! Shared UrbanType instance
  type(UrbanType), public :: urbanxx

end module UrbanxxInstanceMod
