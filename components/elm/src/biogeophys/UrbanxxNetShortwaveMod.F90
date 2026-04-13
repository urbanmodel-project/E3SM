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

end module UrbanxxNetShortwaveMod
