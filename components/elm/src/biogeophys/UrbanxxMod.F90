module UrbanxxMod

  use iso_c_binding
  use urban_mod
  use urban_kokkos_interface
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use spmdMod              , only : masterproc, iam
  use shr_sys_mod          , only : shr_sys_flush
  use elm_varctl           , only : iulog
  use decompMod            , only : bounds_type
  use UrbanParamsType      , only : urbanparams_type
  use SolarAbsorbedType    , only : solarabs_type
  use SurfaceAlbedoType    , only : surfalb_type
  use LandunitType         , only : lun_pp
  use VegetationType       , only : veg_pp
  use FrictionVelocityType , only : frictionvel_type

  implicit none

  private

  type(UrbanType) :: urbanxx

  public :: urbanxx_initialize

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_initialize(bounds, num_urbanl, filter_urbanl, &
       num_urbanc, filter_urbanc, num_urbanp, filter_urbanp, &
       urbanparams_vars, solarabs_vars, surfalb_vars, frictionvel_vars)
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in) :: bounds  
    integer                , intent(in) :: num_urbanl       ! number of urban landunits in clump
    integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
    integer                , intent(in) :: num_urbanc       ! number of urban columns in clump
    integer                , intent(in) :: filter_urbanc(:) ! urban column filter
    integer                , intent(in) :: num_urbanp       ! number of urban patches in clump
    integer                , intent(in) :: filter_urbanp(:) ! urban pft filter
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    type(solarabs_type)    , intent(in) :: solarabs_vars
    type(surfalb_type)     , intent(in) :: surfalb_vars
    type(frictionvel_type) , intent(in) :: frictionvel_vars

    integer :: status

    if (masterproc) then
       write(iulog,*) ' Initializing UrbanXX module...'
    end if

    ! Initialize Kokkos (via C++ wrapper)
    call UrbanKokkosInitialize()

    if (masterproc) then
       write(iulog,*) '=== Fortran Driver with Kokkos ==='
       call UrbanKokkosPrintConfiguration()
    end if

    call UrbanCreate(num_urbanl, urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call SetUrbanParameters(urbanxx, num_urbanl, filter_urbanl, filter_urbanp, &
         urbanparams_vars, frictionvel_vars)

  end subroutine urbanxx_initialize

  !-----------------------------------------------------------------------
  subroutine SetCanyonHwr(urban, num_urbanl, filter_urbanl)
    !
    implicit none
    !
    type(UrbanType) , intent(in)         :: urban
    integer(c_int)  , intent(in)         :: num_urbanl
    integer         , intent(in)         :: filter_urbanl(:) ! urban landunit filter
    !
    integer(c_int)                       :: status
    integer                              ::  fl, l
    real(c_double) , allocatable, target :: canyonHwr(:)

    associate(                           &
         canyon_hwr => lun_pp%canyon_hwr & ! Input:  [real(r8) (:)   ]  ratio of building height to street width
         )
      
      allocate(canyonHwr(num_urbanl))
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         canyonHwr(fl) = lun_pp%canyon_hwr(l)
      end do

      call UrbanSetCanyonHwr(urban, c_loc(canyonHwr), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(*,*) 'Set canyon height-to-width ratio'
      end if

      deallocate(canyonHwr)
    end associate

  end subroutine SetCanyonHwr

  !-----------------------------------------------------------------------
  subroutine SetFracPervRoadOfTotalRoad(urban, num_urbanl, filter_urbanl)
    !
    implicit none
    !
    type(UrbanType) , intent(in)         :: urban
    integer(c_int)  , intent(in)         :: num_urbanl
    integer         , intent(in)         :: filter_urbanl(:) ! urban landunit filter
    !
    integer(c_int)                       :: status
    integer                              :: fl, l
    real(c_double) , allocatable, target :: fracPervRoadOfTotalRoad(:)

    associate(                           &
         wtroad_perv => lun_pp%wtroad_perv & ! Input:  [real(r8) (:)   ]  fraction of pervious road w.r.t. total road
         )

      allocate(fracPervRoadOfTotalRoad(num_urbanl))
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         fracPervRoadOfTotalRoad(fl) = lun_pp%wtroad_perv(l)
      end do

      call UrbanSetFracPervRoadOfTotalRoad(urban, c_loc(fracPervRoadOfTotalRoad), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(*,*) 'Set fraction of pervious road w.r.t. total road'
      end if

      deallocate(fracPervRoadOfTotalRoad)
    end associate

  end subroutine SetFracPervRoadOfTotalRoad

  !-----------------------------------------------------------------------
  subroutine SetWtRoof(urban, num_urbanl, filter_urbanl)
    !
    implicit none
    !
    type(UrbanType) , intent(in)         :: urban
    integer(c_int)  , intent(in)         :: num_urbanl
    integer         , intent(in)         :: filter_urbanl(:) ! urban landunit filter
    !
    integer(c_int)                       :: status
    integer                              :: fl, l
    real(c_double) , allocatable, target :: wtRoof(:)

    associate(                              &
         wtlunit_roof => lun_pp%wtlunit_roof & ! Input:  [real(r8) (:)   ]  roof weight
         )

      allocate(wtRoof(num_urbanl))
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         wtRoof(fl) = lun_pp%wtlunit_roof(l)
      end do

      call UrbanSetWtRoof(urban, c_loc(wtRoof), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(*,*) 'Set roof weight'
      end if

      deallocate(wtRoof)
    end associate

  end subroutine SetWtRoof

  !-----------------------------------------------------------------------
  subroutine SetHeightParameters(urban, num_urbanl, filter_urbanl, filter_urbanp, &
       urbanparams_vars, frictionvel_vars)
    !
    implicit none
    !
    type(UrbanType)        , intent(in) :: urban
    integer(c_int)         , intent(in) :: num_urbanl
    integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
    integer                , intent(in) :: filter_urbanp(:) ! urban pft filter
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    type(frictionvel_type) , intent(in) :: frictionvel_vars
    !
    integer(c_int)                       :: status
    integer                              :: fl, fp, l, p
    real(c_double) , allocatable, target :: forcHgtT(:)
    real(c_double) , allocatable, target :: forcHgtU(:)
    real(c_double) , allocatable, target :: zDTown(:)
    real(c_double) , allocatable, target :: z0Town(:)
    real(c_double) , allocatable, target :: htRoof(:)
    real(c_double) , allocatable, target :: windHgtCanyon(:)

    associate(                                                          &
         forc_hgt_u_patch => frictionvel_vars%forc_hgt_u_patch       , & ! Input: [real(r8) (:)] observational height of wind at pft-level (m)
         forc_hgt_t_patch => frictionvel_vars%forc_hgt_t_patch       , & ! Input: [real(r8) (:)] observational height of temperature at pft-level (m)
         z_d_town         => lun_pp%z_d_town                          , & ! Input: [real(r8) (:)] displacement height of urban landunit (m)
         z_0_town         => lun_pp%z_0_town                          , & ! Input: [real(r8) (:)] momentum roughness length of urban landunit (m)
         ht_roof          => lun_pp%ht_roof                           , & ! Input: [real(r8) (:)] height of urban roof (m)
         wind_hgt_canyon  => urbanparams_vars%wind_hgt_canyon           & ! Input: [real(r8) (:)] height above road at which wind in canyon is to be computed (m)
         )

      allocate(forcHgtT(num_urbanl))
      allocate(forcHgtU(num_urbanl))
      allocate(zDTown(num_urbanl))
      allocate(z0Town(num_urbanl))
      allocate(htRoof(num_urbanl))
      allocate(windHgtCanyon(num_urbanl))

      ! Extract values from patches to landunits - use first urban patch for each landunit
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         forcHgtT(fl) = forc_hgt_t_patch(lun_pp%pfti(l))
         forcHgtU(fl) = forc_hgt_u_patch(lun_pp%pfti(l))
         zDTown(fl) = z_d_town(l)
         z0Town(fl) = z_0_town(l)
         htRoof(fl) = ht_roof(l)
         windHgtCanyon(fl) = wind_hgt_canyon(l)
      end do

      call UrbanSetForcHgtT(urban, c_loc(forcHgtT), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetForcHgtU(urban, c_loc(forcHgtU), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetZDTown(urban, c_loc(zDTown), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetZ0Town(urban, c_loc(z0Town), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetHtRoof(urban, c_loc(htRoof), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetWindHgtCanyon(urban, c_loc(windHgtCanyon), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(*,*) 'Set height parameters'
      end if

      deallocate(forcHgtT)
      deallocate(forcHgtU)
      deallocate(zDTown)
      deallocate(z0Town)
      deallocate(htRoof)
      deallocate(windHgtCanyon)
    end associate

  end subroutine SetHeightParameters

  !-----------------------------------------------------------------------
  subroutine SetUrbanParameters(urban, num_urbanl, filter_urbanl, filter_urbanp, &
       urbanparams_vars, frictionvel_vars)
    !
    implicit none
    !
    type(UrbanType)        , intent(inout) :: urban
    integer(c_int)         , intent(in)    :: num_urbanl
    integer                , intent(in)    :: filter_urbanl(:) ! urban landunit filter
    integer                , intent(in)    :: filter_urbanp(:) ! urban pft filter
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(frictionvel_type) , intent(in)    :: frictionvel_vars

    call SetCanyonHwr(urban, num_urbanl, filter_urbanl)
    call SetFracPervRoadOfTotalRoad(urban, num_urbanl, filter_urbanl)
    call SetWtRoof(urban, num_urbanl, filter_urbanl)
    call SetHeightParameters(urban, num_urbanl, filter_urbanl, filter_urbanp, &
         urbanparams_vars, frictionvel_vars)

  end subroutine SetUrbanParameters

end module UrbanxxMod
