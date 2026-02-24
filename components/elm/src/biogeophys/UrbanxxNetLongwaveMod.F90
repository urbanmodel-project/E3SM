module UrbanxxNetLongwaveMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute net longwave radiation for the Urban++ model.
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
  use ColumnDataType       , only : col_es, col_pp
  use UrbanxxInstanceMod   , only : urbanxx

  implicit none

  private

  ! Persistent input buffers (allocated once in init)
  real(c_double) , allocatable, target :: t_roof(:)
  real(c_double) , allocatable, target :: t_improad(:)
  real(c_double) , allocatable, target :: t_perroad(:)
  real(c_double) , allocatable, target :: t_sunwall(:)
  real(c_double) , allocatable, target :: t_shadwall(:)

  ! Persistent output buffers (allocated once in init)
  real(c_double) , allocatable, target, public :: lwnet_roof(:)
  real(c_double) , allocatable, target, public :: lwnet_improad(:)
  real(c_double) , allocatable, target, public :: lwnet_perroad(:)
  real(c_double) , allocatable, target, public :: lwnet_sunwall(:)
  real(c_double) , allocatable, target, public :: lwnet_shadwall(:)
  real(c_double) , allocatable, target, public :: lwup_roof(:)
  real(c_double) , allocatable, target, public :: lwup_improad(:)
  real(c_double) , allocatable, target, public :: lwup_perroad(:)
  real(c_double) , allocatable, target, public :: lwup_sunwall(:)
  real(c_double) , allocatable, target, public :: lwup_shadwall(:)

  public :: urbanxx_netLongwave_init
  public :: urbanxx_netLongwave

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_netLongwave_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for net longwave computation.
    ! Called once during initialization.
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl

    ! Input buffers
    allocate(t_roof(num_urbanl))
    allocate(t_improad(num_urbanl))
    allocate(t_perroad(num_urbanl))
    allocate(t_sunwall(num_urbanl))
    allocate(t_shadwall(num_urbanl))

    ! Output buffers
    allocate(lwnet_roof(num_urbanl))
    allocate(lwnet_improad(num_urbanl))
    allocate(lwnet_perroad(num_urbanl))
    allocate(lwnet_sunwall(num_urbanl))
    allocate(lwnet_shadwall(num_urbanl))
    allocate(lwup_roof(num_urbanl))
    allocate(lwup_improad(num_urbanl))
    allocate(lwup_perroad(num_urbanl))
    allocate(lwup_sunwall(num_urbanl))
    allocate(lwup_shadwall(num_urbanl))

  end subroutine urbanxx_netLongwave_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_netLongwave(num_urbanl, filter_urbanl, surfalb_vars, &
         urbanparams_vars, frictionvel_vars)
    !
    use column_varcon       , only : icol_road_perv, icol_road_imperv
    use column_varcon       , only : icol_roof, icol_sunwall, icol_shadewall
    !
    implicit none
    !
    integer(c_int)         , intent(in) :: num_urbanl
    integer                , intent(in) :: filter_urbanl(:)         ! urban landunit filter
    type(surfalb_type)     , intent(in) :: surfalb_vars
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    type(frictionvel_type) , intent(in) :: frictionvel_vars
    !
    integer                              :: fl, l, c_start, c_end, c
    integer(c_int)                       :: status

    associate(                       &
         ctype  =>    col_pp%itype , & ! Input:  [integer (:)    ]  column type
         coli   =>    lun_pp%coli  , & ! Input:  [integer (:)    ]  beginning column index for landunit
         colf   =>    lun_pp%colf  , & ! Input:  [integer (:)    ]  ending column index for landunit
         t_grnd =>    col_es%t_grnd  & ! Input:  [real(r8) (:)   ]  ground temperature (K)
         )

      ! Extract surface temperatures from columns to landunits
      ! For urban landunits, there are multiple columns per landunit representing different surfaces.
      ! We extract temperatures based on column types.
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         c_start = coli(l)
         c_end   = colf(l)

         do c = c_start, c_end
            select case (ctype(c))
            case (icol_roof)  ! Roof
               t_roof(fl) = t_grnd(c)
            case (icol_road_imperv)  ! Impervious Road
               t_improad(fl) = t_grnd(c)
            case (icol_road_perv)  ! Pervious Road
               t_perroad(fl) = t_grnd(c)
            case (icol_sunwall)  ! Sunlit Wall
               t_sunwall(fl) = t_grnd(c)
            case (icol_shadewall)  ! Shaded Wall
               t_shadwall(fl) = t_grnd(c)
            case default
               ! Do nothing for other types
            end select
         end do
      end do

      ! set surface temperatures in UrbanXX
      call UrbanSetEffectiveSurfTempRoof(urbanxx, c_loc(t_roof), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetEffectiveSurfTempImperviousRoad(urbanxx, c_loc(t_improad), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetEffectiveSurfTempPerviousRoad(urbanxx, c_loc(t_perroad), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetEffectiveSurfTempSunlitWall(urbanxx, c_loc(t_sunwall), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetEffectiveSurfTempShadedWall(urbanxx, c_loc(t_shadwall), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanComputeNetLongwave(urbanxx, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! Extract net longwave radiation from UrbanXX
      call UrbanGetNetLongwaveRoof(urbanxx, c_loc(lwnet_roof), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetNetLongwaveImperviousRoad(urbanxx, c_loc(lwnet_improad), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetNetLongwavePerviousRoad(urbanxx, c_loc(lwnet_perroad), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetNetLongwaveSunlitWall(urbanxx, c_loc(lwnet_sunwall), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetNetLongwaveShadedWall(urbanxx, c_loc(lwnet_shadwall), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! Extract upward longwave radiation from UrbanXX
      call UrbanGetUpwardLongwaveRoof(urbanxx, c_loc(lwup_roof), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetUpwardLongwaveImperviousRoad(urbanxx, c_loc(lwup_improad), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetUpwardLongwavePerviousRoad(urbanxx, c_loc(lwup_perroad), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetUpwardLongwaveSunlitWall(urbanxx, c_loc(lwup_sunwall), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetUpwardLongwaveShadedWall(urbanxx, c_loc(lwup_shadwall), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    end associate

  end subroutine urbanxx_netLongwave

end module UrbanxxNetLongwaveMod
