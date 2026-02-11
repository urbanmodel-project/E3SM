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
  use TopounitDataType     , only : top_as, top_af
  use ColumnDataType       , only : col_es, col_pp
  use LandunitDataType     , only : lun_es, lun_ws
  use SoilStateType        , only : soilstate_type
  use abortutils           , only : endrun

  implicit none

  private

  ! Constants for radiation bands and types
  integer(c_int), parameter :: numBands = 2  ! VIS, NIR
  integer(c_int), parameter :: numTypes = 2  ! Direct, Diffuse

  type(UrbanType) :: urbanxx

  public :: urbanxx_initialize
  public :: urbanxx_SetAtmosphericForcing
  public :: urbanxx_netLongwave
  public :: urbanxx_netShortwave
  public :: urbanxx_surfaceFluxes
contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_initialize(bounds, num_urbanl, filter_urbanl, &
       num_urbanc, filter_urbanc, num_urbanp, filter_urbanp, &
       urbanparams_vars, solarabs_vars, surfalb_vars, frictionvel_vars, &
       soilstate_vars)
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
    type(soilstate_type)   , intent(in) :: soilstate_vars

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

    call SetUrbanParameters(urbanxx, num_urbanl, filter_urbanl, &
         urbanparams_vars, frictionvel_vars, soilstate_vars)

    ! Setup urban model (initialize temperatures and other setup tasks)
    call UrbanSetup(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    if (masterproc) then
       write(*,*) 'Completed urban model setup'
    end if

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
  subroutine SetHeightParameters(urban, num_urbanl, filter_urbanl, &
       urbanparams_vars, frictionvel_vars)
    !
    implicit none
    !
    type(UrbanType)        , intent(in) :: urban
    integer(c_int)         , intent(in) :: num_urbanl
    integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
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
         forcHgtT(fl)      = forc_hgt_t_patch(lun_pp%pfti(l))
         forcHgtU(fl)      = forc_hgt_u_patch(lun_pp%pfti(l))
         zDTown(fl)        = z_d_town(l)
         z0Town(fl)        = z_0_town(l)
         htRoof(fl)        = ht_roof(l)
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
  subroutine SetAlbedo(urban, num_urbanl, filter_urbanl, urbanparams_vars)
    !
    implicit none
    !
    type(UrbanType)        , intent(in) :: urban
    integer(c_int)         , intent(in) :: num_urbanl
    integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    !
    integer(c_int)                       :: status
    integer                              :: fl, l, iband, itype, count
    integer(c_int)                       :: totalSize3D
    integer(c_int), dimension(3)         :: size3D
    real(c_double) , allocatable, target :: albedoPerviousRoad(:)
    real(c_double) , allocatable, target :: albedoImperviousRoad(:)
    real(c_double) , allocatable, target :: albedoSunlitWall(:)
    real(c_double) , allocatable, target :: albedoShadedWall(:)
    real(c_double) , allocatable, target :: albedoRoof(:)

    associate(                                                    &
         alb_roof_dir       => urbanparams_vars%alb_roof_dir    , & ! Input: [real(r8) (:,:)] direct roof albedo
         alb_roof_dif       => urbanparams_vars%alb_roof_dif    , & ! Input: [real(r8) (:,:)] diffuse roof albedo
         alb_improad_dir    => urbanparams_vars%alb_improad_dir , & ! Input: [real(r8) (:,:)] direct impervious road albedo
         alb_improad_dif    => urbanparams_vars%alb_improad_dif , & ! Input: [real(r8) (:,:)] diffuse imprevious road albedo
         alb_perroad_dir    => urbanparams_vars%alb_perroad_dir , & ! Input: [real(r8) (:,:)] direct pervious road albedo
         alb_perroad_dif    => urbanparams_vars%alb_perroad_dif , & ! Input: [real(r8) (:,:)] diffuse pervious road albedo
         alb_wall_dir       => urbanparams_vars%alb_wall_dir    , & ! Input: [real(r8) (:,:)] direct wall albedo
         alb_wall_dif       => urbanparams_vars%alb_wall_dif      & ! Input: [real(r8) (:,:)] diffuse wall albedo
         )

      size3D = [num_urbanl, numBands, numTypes]
      totalSize3D = num_urbanl * numBands * numTypes

      allocate(albedoPerviousRoad(totalSize3D))
      allocate(albedoImperviousRoad(totalSize3D))
      allocate(albedoSunlitWall(totalSize3D))
      allocate(albedoShadedWall(totalSize3D))
      allocate(albedoRoof(totalSize3D))

      ! Fill arrays using same indexing as C: idx = ilandunit * numBands * numTypes + iband * numTypes + itype
      ! itype = 0 corresponds to diffuse (*_dif), itype = 1 corresponds to direct (*_dir)
      ! Note: Fortran arrays are 1-indexed, so we adjust accordingly
      !
      ! In Urbanxx lib, the albedos are defined as following 3D: albedo(l, ib, itype)
      
      if (UrbanKokkosIsLayoutLeft()) then

         ! Kokkos layout is left: Thus, the 'l'-index will be incremented fastest

         count = 0
         do itype = 0, 1
            do iband = 0, numBands - 1
               do fl = 1, num_urbanl
                  l = filter_urbanl(fl)
                  count = count + 1
                  if (itype == 0) then
                     ! itype = 0: diffuse
                     albedoPerviousRoad(count)   = alb_perroad_dif(l, iband+1)
                     albedoImperviousRoad(count) = alb_improad_dif(l, iband+1)
                     albedoSunlitWall(count)     = alb_wall_dif(l, iband+1)
                     albedoShadedWall(count)     = alb_wall_dif(l, iband+1)
                     albedoRoof(count)           = alb_roof_dif(l, iband+1)

                  else
                     ! itype = 1: direct
                     albedoPerviousRoad(count)   = alb_perroad_dir(l, iband+1)
                     albedoImperviousRoad(count) = alb_improad_dir(l, iband+1)
                     albedoSunlitWall(count)     = alb_wall_dir(l, iband+1)
                     albedoShadedWall(count)     = alb_wall_dir(l, iband+1)
                     albedoRoof(count)           = alb_roof_dir(l, iband+1)
                  endif
               enddo
            enddo
         enddo
      else

         ! Kokkos layout is right: Thus, the 'itype'-index will be incremented fastest

         count = 0
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            do iband = 0, numBands - 1
               do itype = 0, 1
                  count = count + 1
                  if (itype == 0) then
                     ! itype = 0: diffuse
                     albedoPerviousRoad(count)   = alb_perroad_dif(l, iband+1)
                     albedoImperviousRoad(count) = alb_improad_dif(l, iband+1)
                     albedoSunlitWall(count)     = alb_wall_dif(l, iband+1)
                     albedoShadedWall(count)     = alb_wall_dif(l, iband+1)
                     albedoRoof(count)           = alb_roof_dif(l, iband+1)
                  else
                     ! itype = 1: direct
                     albedoPerviousRoad(count)   = alb_perroad_dir(l, iband+1)
                     albedoImperviousRoad(count) = alb_improad_dir(l, iband+1)
                     albedoSunlitWall(count)     = alb_wall_dir(l, iband+1)
                     albedoShadedWall(count)     = alb_wall_dir(l, iband+1)
                     albedoRoof(count)           = alb_roof_dir(l, iband+1)
                  endif
               enddo
            enddo
         enddo
      end if

      call UrbanSetAlbedoPerviousRoad(urban, c_loc(albedoPerviousRoad), size3D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAlbedoImperviousRoad(urban, c_loc(albedoImperviousRoad), size3D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAlbedoSunlitWall(urban, c_loc(albedoSunlitWall), size3D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAlbedoShadedWall(urban, c_loc(albedoShadedWall), size3D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAlbedoRoof(urban, c_loc(albedoRoof), size3D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(*,*) 'Set albedo values for all surfaces'
      end if

      deallocate(albedoPerviousRoad)
      deallocate(albedoImperviousRoad)
      deallocate(albedoSunlitWall)
      deallocate(albedoShadedWall)
      deallocate(albedoRoof)
    end associate

  end subroutine SetAlbedo

  !-----------------------------------------------------------------------
  subroutine SetEmissivity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
    !
    implicit none
    !
    type(UrbanType)        , intent(in) :: urban
    integer(c_int)         , intent(in) :: num_urbanl
    integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    !
    integer(c_int)                       :: status
    integer                              :: fl, l
    real(c_double) , allocatable, target :: emissivityPerviousRoad(:)
    real(c_double) , allocatable, target :: emissivityImperviousRoad(:)
    real(c_double) , allocatable, target :: emissivityWall(:)
    real(c_double) , allocatable, target :: emissivityRoof(:)

    associate(                                       &
         em_roof    => urbanparams_vars%em_roof    , & ! Input: [real(r8) (:)] roof emissivity
         em_improad => urbanparams_vars%em_improad , & ! Input: [real(r8) (:)] impervious road emissivity
         em_perroad => urbanparams_vars%em_perroad , & ! Input: [real(r8) (:)] pervious road emissivity
         em_wall    => urbanparams_vars%em_wall      & ! Input: [real(r8) (:)] wall emissivity
         )

      allocate(emissivityPerviousRoad(num_urbanl))
      allocate(emissivityImperviousRoad(num_urbanl))
      allocate(emissivityWall(num_urbanl))
      allocate(emissivityRoof(num_urbanl))

      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         emissivityPerviousRoad(fl)   = em_perroad(l)
         emissivityImperviousRoad(fl) = em_improad(l)
         emissivityWall(fl)           = em_wall(l)
         emissivityRoof(fl)           = em_roof(l)
      end do

      call UrbanSetEmissivityPerviousRoad(urban, c_loc(emissivityPerviousRoad), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetEmissivityImperviousRoad(urban, c_loc(emissivityImperviousRoad), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetEmissivityWall(urban, c_loc(emissivityWall), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetEmissivityRoof(urban, c_loc(emissivityRoof), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(*,*) 'Set emissivity values for all surfaces'
      end if

      deallocate(emissivityPerviousRoad)
      deallocate(emissivityImperviousRoad)
      deallocate(emissivityWall)
      deallocate(emissivityRoof)
    end associate

  end subroutine SetEmissivity

  !-----------------------------------------------------------------------
  subroutine SetThermalConductivity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
    !
    use elm_varpar, only : nlevgrnd, nlevurb
    !
    implicit none
    !
    type(UrbanType)        , intent(in) :: urban
    integer(c_int)         , intent(in) :: num_urbanl
    integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    !
    integer(c_int)                       :: status
    integer                              :: fl, l, j, idx
    real(c_double) , allocatable, target :: tkRoad(:)
    real(c_double) , allocatable, target :: tkWall(:)
    real(c_double) , allocatable, target :: tkRoof(:)
    integer(c_int), dimension(2) :: size2D_road, size2D_urban

    associate(                                       &
         tk_wall    => urbanparams_vars%tk_wall    , & ! Input: [real(r8) (:,:)] thermal conductivity of urban wall
         tk_roof    => urbanparams_vars%tk_roof    , & ! Input: [real(r8) (:,:)] thermal conductivity of urban roof
         tk_improad => urbanparams_vars%tk_improad   & ! Input: [real(r8) (:,:)] thermal conductivity of urban impervious road
         )

      size2D_road(1) = num_urbanl
      size2D_road(2) = nlevgrnd

      size2D_urban(1) = num_urbanl
      size2D_urban(2) = nlevurb
      allocate(tkRoad(num_urbanl * nlevgrnd))
      allocate(tkWall(num_urbanl * nlevurb))
      allocate(tkRoof(num_urbanl * nlevurb))

      ! URBANXX_FIX_ME: Currently using only first layer values.
      ! ELM has multi-layer thermal conductivity data (tk_wall, tk_roof, tk_improad are dimensioned as [landunit, nlevurb]).
      ! Urban++ may need to be updated to accept multi-layer thermal properties.
      idx = 0
      do j = 1, nlevgrnd
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            idx = idx + 1
            tkRoad(idx) = tk_improad(l, j)
         enddo
      end do

      idx = 0
      do j = 1, num_urbanl
         do fl = 1, nlevurb
            l = filter_urbanl(fl)
            idx = idx + 1
            tkWall(idx) = tk_wall(l, j)
            tkRoof(idx) = tk_roof(l, j)
         enddo
      end do

      call UrbanSetThermalConductivityRoad(urban, c_loc(tkRoad), &
           size2D_road, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetThermalConductivityWall(urban, c_loc(tkWall), &
           size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetThermalConductivityRoof(urban, c_loc(tkRoof), &
           size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(iulog,*) 'Set thermal conductivity values for all surfaces'
         write(iulog,*) 'URBANXX_FIX_ME: Only using first layer of multi-layer thermal conductivity data'
      end if

      deallocate(tkRoad)
      deallocate(tkWall)
      deallocate(tkRoof)
    end associate

  end subroutine SetThermalConductivity

  !-----------------------------------------------------------------------
  subroutine SetHeatCapacity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
    !
    use elm_varpar, only : nlevgrnd, nlevurb
    !
    implicit none
    !
    type(UrbanType)        , intent(in) :: urban
    integer(c_int)         , intent(in) :: num_urbanl
    integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    !
    integer(c_int)                       :: status
    integer                              :: fl, l, j, idx
    real(c_double) , allocatable, target :: cvRoad(:)
    real(c_double) , allocatable, target :: cvWall(:)
    real(c_double) , allocatable, target :: cvRoof(:)
    integer(c_int), dimension(2) :: size2D_road, size2D_urban

    associate(                                           &
         cv_wall      => urbanparams_vars%cv_wall      , & ! Input: [real(r8) (:,:)] heat capacity of urban wall
         cv_roof      => urbanparams_vars%cv_roof      , & ! Input: [real(r8) (:,:)] heat capacity of urban roof
         cv_improad   => urbanparams_vars%cv_improad     & ! Input: [real(r8) (:,:)] heat capacity of urban impervious road
         )

      size2D_road(1) = num_urbanl
      size2D_road(2) = nlevgrnd

      size2D_urban(1) = num_urbanl
      size2D_urban(2) = nlevurb

      allocate(cvRoad(num_urbanl * nlevgrnd))
      allocate(cvWall(num_urbanl * nlevurb))
      allocate(cvRoof(num_urbanl * nlevurb))

      idx = 0
      do j = 1, nlevgrnd
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            idx = idx + 1
            cvRoad(idx) = cv_improad(l, j)
         enddo
      end do

      idx = 0
      do j = 1, nlevurb
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            idx = idx + 1
            cvWall(idx) = cv_wall(l, j)
            cvRoof(idx) = cv_roof(l, j)
         enddo
      end do

      call UrbanSetHeatCapacityRoad(urban, c_loc(cvRoad), &
           size2D_road, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetHeatCapacityWall(urban, c_loc(cvWall), &
           size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetHeatCapacityRoof(urban, c_loc(cvRoof), &
           size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(iulog,*) 'Set heat capacity values for all surfaces'
      end if

      deallocate(cvRoad)
      deallocate(cvWall)
      deallocate(cvRoof)
    end associate

  end subroutine SetHeatCapacity

  !-----------------------------------------------------------------------
  subroutine urbanxx_SetAtmosphericForcing(num_urbanl, filter_urbanl, surfalb_vars, &
         urbanparams_vars, frictionvel_vars)
    !
    implicit none
    !
    integer(c_int)     , intent(in) :: num_urbanl
    integer            , intent(in) :: filter_urbanl(:)         ! urban landunit filter
    type(surfalb_type) , intent(in) :: surfalb_vars
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(frictionvel_type) , intent(in)    :: frictionvel_vars

    !
    integer(c_int)                       :: status
    integer                              :: fl, l, t, iband, itype, idx
    integer(c_int)                       :: totalSize3D
    integer(c_int), dimension(3)         :: size3D
    real(c_double) , allocatable, target :: atmTemp(:)
    real(c_double) , allocatable, target :: atmPotTemp(:)
    real(c_double) , allocatable, target :: atmRho(:)
    real(c_double) , allocatable, target :: atmSpcHumd(:)
    real(c_double) , allocatable, target :: atmPress(:)
    real(c_double) , allocatable, target :: atmWindU(:)
    real(c_double) , allocatable, target :: atmWindV(:)
    real(c_double) , allocatable, target :: atmCoszen(:)
    real(c_double) , allocatable, target :: atmFracSnow(:)
    real(c_double) , allocatable, target :: atmLongwave(:)
    real(c_double) , allocatable, target :: atmShortwave(:)

    associate(                                                        &
         forc_t     => top_as%tbot     , & ! Input: [real(r8) (:)] atmospheric temperature (K)
         forc_th    => top_as%thbot    , & ! Input: [real(r8) (:)] atmospheric potential temperature (K)
         forc_rho   => top_as%rhobot   , & ! Input: [real(r8) (:)] air density (kg/m**3)
         forc_q     => top_as%qbot     , & ! Input: [real(r8) (:)] atmospheric specific humidity (kg/kg)
         forc_pbot  => top_as%pbot     , & ! Input: [real(r8) (:)] atmospheric pressure (Pa)
         forc_u     => top_as%ubot     , & ! Input: [real(r8) (:)] atmospheric wind speed in east direction (m/s)
         forc_v     => top_as%vbot     , & ! Input: [real(r8) (:)] atmospheric wind speed in north direction (m/s)
         forc_lwrad => top_af%lwrad_pp , & ! Input: [real(r8) (:)] downward infrared (longwave) radiation under PP (W/m**2)
         forc_snow  => top_af%snow     , & ! Input: [real(r8) (:)] downscaled snow
         forc_solad => top_af%solad_pp , & ! Input: [real(r8) (:,:)] direct beam radiation under PP (vis=forc_sols , nir=forc_soll ) (W/m**2)
         forc_solai => top_af%solai_pp , & ! Input: [real(r8) (:,:)] diffuse beam radiation under PP (vis=forc_sols , nir=forc_soll ) (W/m**2)
         coli       => lun_pp%coli                           & ! Input: [integer (:)] beginning column index for landunit
         )

      ! Allocate arrays
      allocate(atmTemp(num_urbanl))
      allocate(atmPotTemp(num_urbanl))
      allocate(atmRho(num_urbanl))
      allocate(atmSpcHumd(num_urbanl))
      allocate(atmPress(num_urbanl))
      allocate(atmWindU(num_urbanl))
      allocate(atmWindV(num_urbanl))
      allocate(atmCoszen(num_urbanl))
      allocate(atmFracSnow(num_urbanl))
      allocate(atmLongwave(num_urbanl))

      size3D = [num_urbanl, numBands, numTypes]
      totalSize3D = num_urbanl * numBands * numTypes
      allocate(atmShortwave(totalSize3D))

      ! Fill arrays with values from ELM data structures
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         t = lun_pp%topounit(l)

         atmTemp(fl)     = forc_t(t)
         atmPotTemp(fl)  = forc_th(t)
         atmRho(fl)      = forc_rho(t)
         atmSpcHumd(fl)  = forc_q(t)
         atmPress(fl)    = forc_pbot(t)
         atmWindU(fl)    = forc_u(t)
         atmWindV(fl)    = forc_v(t)
         atmCoszen(fl)   = surfalb_vars%coszen_col(coli(l))  ! Assumes coszen for each column are the same
         atmFracSnow(fl) = forc_snow(t)
         atmLongwave(fl) = forc_lwrad(t)
      end do

      ! Fill shortwave arrays with direct and diffuse for VIS and NIR bands
      ! Indexing: idx = ilandunit * numBands * numTypes + iband * numTypes + itype
      ! itype = 0: diffuse, itype = 1: direct
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         t = lun_pp%topounit(l)
         do iband = 0, numBands - 1
            ! itype = 0: diffuse
            idx = (fl-1) * numBands * numTypes + iband * numTypes + 0 + 1  ! +1 for Fortran 1-indexing
            atmShortwave(idx) = forc_solai(t, iband+1)

            ! itype = 1: direct
            idx = (fl-1) * numBands * numTypes + iband * numTypes + 1 + 1  ! +1 for Fortran 1-indexing
            atmShortwave(idx) = forc_solad(t, iband+1)
         end do
      end do

      ! Set atmospheric forcing
      call UrbanSetAtmTemp(urbanxx, c_loc(atmTemp), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmPotTemp(urbanxx, c_loc(atmPotTemp), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmRho(urbanxx, c_loc(atmRho), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmSpcHumd(urbanxx, c_loc(atmSpcHumd), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmPress(urbanxx, c_loc(atmPress), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmWindU(urbanxx, c_loc(atmWindU), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmWindV(urbanxx, c_loc(atmWindV), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmCoszen(urbanxx, c_loc(atmCoszen), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmFracSnow(urbanxx, c_loc(atmFracSnow), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmLongwaveDown(urbanxx, c_loc(atmLongwave), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanSetAtmShortwaveDown(urbanxx, c_loc(atmShortwave), size3D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(iulog,*) 'Set atmospheric forcing from ELM data structures'
      end if

      call SetHeightParameters(urbanxx, num_urbanl, filter_urbanl, &
         urbanparams_vars, frictionvel_vars)

      ! Free arrays
      deallocate(atmTemp)
      deallocate(atmPotTemp)
      deallocate(atmRho)
      deallocate(atmSpcHumd)
      deallocate(atmPress)
      deallocate(atmWindU)
      deallocate(atmWindV)
      deallocate(atmCoszen)
      deallocate(atmFracSnow)
      deallocate(atmLongwave)
      deallocate(atmShortwave)
    end associate

  end subroutine urbanxx_SetAtmosphericForcing

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
    real(c_double) , allocatable, target :: t_roof(:)
    real(c_double) , allocatable, target :: t_improad(:)
    real(c_double) , allocatable, target :: t_perroad(:)
    real(c_double) , allocatable, target :: t_sunwall(:)
    real(c_double) , allocatable, target :: t_shadwall(:)

    associate(                       &
         ctype  =>    col_pp%itype , & ! Input:  [integer (:)    ]  column type
         coli   =>    lun_pp%coli  , & ! Input:  [integer (:)    ]  beginning column index for landunit
         colf   =>    lun_pp%colf  , & ! Input:  [integer (:)    ]  ending column index for landunit
         t_grnd =>    col_es%t_grnd  & ! Input:  [real(r8) (:)   ]  ground temperature (K)
         )

      allocate(t_roof(num_urbanl))
      allocate(t_improad(num_urbanl))
      allocate(t_perroad(num_urbanl))
      allocate(t_sunwall(num_urbanl))
      allocate(t_shadwall(num_urbanl))

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

      ! Free arrays
      deallocate(t_roof)
      deallocate(t_improad)
      deallocate(t_perroad)
      deallocate(t_sunwall)
      deallocate(t_shadwall)

    end associate

   end subroutine urbanxx_netLongwave

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

   !-----------------------------------------------------------------------
   subroutine urbanxx_surfaceFluxes(num_urbanl, filter_urbanl, surfalb_vars, &
        urbanparams_vars, frictionvel_vars)
     !
     implicit none
     !
     integer(c_int)     , intent(in) :: num_urbanl
     integer            , intent(in) :: filter_urbanl(:)         ! urban landunit filter
     type(surfalb_type) , intent(in) :: surfalb_vars
     type(urbanparams_type) , intent(in)    :: urbanparams_vars
     type(frictionvel_type) , intent(in)    :: frictionvel_vars
     !
     integer(c_int)                       :: status

     call SetHeightParameters(urbanxx, num_urbanl, filter_urbanl, &
        urbanparams_vars, frictionvel_vars)
     call UrbanComputeSurfaceFluxes(urbanxx, status)
     if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

   end subroutine urbanxx_surfaceFluxes

   !-----------------------------------------------------------------------
   subroutine SetSoilProperties(urban, num_urbanl, filter_urbanl, soilstate_vars)
     !
     use elm_varpar, only : nlevgrnd
     !
     implicit none
     !
     type(UrbanType)      , intent(in)    :: urban
     integer(c_int)       , intent(in)    :: num_urbanl
     integer              , intent(in)    :: filter_urbanl(:) ! urban landunit filter
     type(soilstate_type) , intent(in)    :: soilstate_vars
     !
     integer(c_int)                       :: status
     integer                              :: fl, l, c, j, idx
     integer(c_int)                       :: totalSize
     integer(c_int), dimension(2)         :: size2D
     logical(c_bool)                      :: isLayoutLeft
     real(c_double) , allocatable, target :: sand(:)
     real(c_double) , allocatable, target :: clay(:)
     real(c_double) , allocatable, target :: organic(:)

     associate(                               &
          cellsand => soilstate_vars%cellsand_col , & ! Input: [real(r8) (:,:)] sand fraction
          cellclay => soilstate_vars%cellclay_col , & ! Input: [real(r8) (:,:)] clay fraction
          cellorg  => soilstate_vars%cellorg_col  , & ! Input: [real(r8) (:,:)] organic matter
          coli     => lun_pp%coli                   & ! Input: [integer (:)] beginning column index for landunit
          )

       totalSize = num_urbanl * nlevgrnd
       size2D(1) = num_urbanl
       size2D(2) = nlevgrnd

       allocate(sand(totalSize))
       allocate(clay(totalSize))
       allocate(organic(totalSize))

       ! Check Kokkos memory layout
       isLayoutLeft = UrbanKokkosIsLayoutLeft()

       if (isLayoutLeft) then
         ! LayoutLeft: First dimension (landunits) varies fastest
         ! Iterate: layer (outer), landunits (inner)
         idx = 0
         do j = 1, nlevgrnd
           do fl = 1, num_urbanl
             l = filter_urbanl(fl)
             c = coli(l)  ! Get first column for this landunit
             idx = idx + 1
             sand(idx) = cellsand(c, j)
             clay(idx) = cellclay(c, j)
             organic(idx) = cellorg(c, j)
           end do
         end do
       else
         ! LayoutRight: Last dimension (layers) varies fastest
         ! Iterate: landunits (outer), layer (inner)
         idx = 0
         do fl = 1, num_urbanl
           l = filter_urbanl(fl)
           c = coli(l)  ! Get first column for this landunit
           do j = 1, nlevgrnd
             idx = idx + 1
             sand(idx) = cellsand(c, j)
             clay(idx) = cellclay(c, j)
             organic(idx) = cellorg(c, j)
           end do
         end do
       end if

       call UrbanSetSandPerviousRoad(urban, c_loc(sand), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetClayPerviousRoad(urban, c_loc(clay), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetOrganicPerviousRoad(urban, c_loc(organic), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       if (masterproc) then
         write(iulog,*) 'Set soil properties for pervious road (sand, clay, organic) from ELM data'
       end if

       deallocate(sand)
       deallocate(clay)
       deallocate(organic)
     end associate

   end subroutine SetSoilProperties

   !-----------------------------------------------------------------------
   subroutine SetCanyonAirStates(urban, num_urbanl, filter_urbanl)
     !
     ! !DESCRIPTION:
     ! Set canyon air temperature and specific humidity from landunit state data.
     !
     implicit none
     !
     type(UrbanType), intent(in) :: urban
     integer(c_int) , intent(in) :: num_urbanl
     integer        , intent(in) :: filter_urbanl(:) ! urban landunit filter
     !
     integer(c_int)                       :: status
     integer                              :: fl, l
     real(c_double) , allocatable, target :: tempCanyonAir(:)
     real(c_double) , allocatable, target :: qafCanyonAir(:)

     associate(                     &
          taf => lun_es%taf      , & ! Input: [real(r8) (:)] urban canopy air temperature (K)
          qaf => lun_ws%qaf        & ! Input: [real(r8) (:)] urban canopy air specific humidity (kg H2O/kg moist air)
          )

       ! Allocate arrays
       allocate(tempCanyonAir(num_urbanl))
       allocate(qafCanyonAir(num_urbanl))

       ! Copy data from landunit state to urban arrays
       do fl = 1, num_urbanl
          l = filter_urbanl(fl)
          tempCanyonAir(fl) = taf(l)
          qafCanyonAir(fl) = qaf(l)
       end do

       ! Set canyon air temperature
       call UrbanSetCanyonAirTemperature(urban, c_loc(tempCanyonAir), num_urbanl, status)
       if (status /= URBAN_SUCCESS) then
          write(iulog,*) 'ERROR: UrbanSetCanyonAirTemperature failed with status ', status
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       ! Set canyon specific humidity
       call UrbanSetCanyonSpecificHumidity(urban, c_loc(qafCanyonAir), num_urbanl, status)
       if (status /= URBAN_SUCCESS) then
          write(iulog,*) 'ERROR: UrbanSetCanyonSpecificHumidity failed with status ', status
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       ! Deallocate arrays
       deallocate(tempCanyonAir)
       deallocate(qafCanyonAir)

     end associate

   end subroutine SetCanyonAirStates

   !-----------------------------------------------------------------------
   subroutine SetNumberOfActiveLayersImperviousRoad(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     !
     implicit none
     !
     type(UrbanType)        , intent(in) :: urban
     integer(c_int)         , intent(in) :: num_urbanl
     integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
     type(urbanparams_type) , intent(in) :: urbanparams_vars
     !
     integer(c_int)                       :: status
     integer                              :: fl, l
     real(c_double) , allocatable, target :: numActiveLayers(:)

     associate(                                 &
          nlev_improad => urbanparams_vars%nlev_improad  & ! Input: [integer (:)] number of impervious road layers
          )

       ! Allocate and populate the number of active layers array
       allocate(numActiveLayers(num_urbanl))
       do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         numActiveLayers(fl) = real(nlev_improad(l), c_double)
       end do

       ! Set the number of active layers in the Urban++ instance
       call UrbanSetNumberOfActiveLayersImperviousRoad(urban, &
            c_loc(numActiveLayers), num_urbanl, status)
       if (status /= URBAN_SUCCESS) then
         write(iulog,*) 'ERROR: UrbanSetNumberOfActiveLayersImperviousRoad failed with status: ', status
         call UrbanError(iam, __LINE__, status)
       end if

       deallocate(numActiveLayers)

     end associate

   end subroutine SetNumberOfActiveLayersImperviousRoad

   !-----------------------------------------------------------------------
   subroutine SetBuildingTemperature(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     !
     implicit none
     !
     type(UrbanType)        , intent(in) :: urban
     integer(c_int)         , intent(in) :: num_urbanl
     integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
     type(urbanparams_type) , intent(in) :: urbanparams_vars
     !
     integer(c_int)                       :: status
     integer                              :: fl, l
     real(c_double) , allocatable, target :: minTemp(:)
     real(c_double) , allocatable, target :: maxTemp(:)
     real(c_double) , allocatable, target :: wallThickness(:)
     real(c_double) , allocatable, target :: roofThickness(:)

     associate(                                        &
          t_building_min => urbanparams_vars%t_building_min , & ! Input: [real(r8) (:)] minimum internal building temperature (K)
          t_building_max => urbanparams_vars%t_building_max , & ! Input: [real(r8) (:)] maximum internal building temperature (K)
          thick_wall     => urbanparams_vars%thick_wall     , & ! Input: [real(r8) (:)] total thickness of urban wall (m)
          thick_roof     => urbanparams_vars%thick_roof       & ! Input: [real(r8) (:)] total thickness of urban roof (m)
          )

       ! Set building temperature limits
       allocate(minTemp(num_urbanl))
       allocate(maxTemp(num_urbanl))

       do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         minTemp(fl) = t_building_min(l)
         maxTemp(fl) = t_building_max(l)
       end do

       call UrbanSetBuildingMinTemperature(urban, c_loc(minTemp), &
            num_urbanl, status)
       if (status /= URBAN_SUCCESS) then
         write(iulog,*) 'ERROR: UrbanSetBuildingMinTemperature failed with status: ', status
         call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       call UrbanSetBuildingMaxTemperature(urban, c_loc(maxTemp), &
            num_urbanl, status)
       if (status /= URBAN_SUCCESS) then
         write(iulog,*) 'ERROR: UrbanSetBuildingMaxTemperature failed with status: ', status
         call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       deallocate(minTemp)
       deallocate(maxTemp)

       ! Set building thickness parameters
       allocate(wallThickness(num_urbanl))
       allocate(roofThickness(num_urbanl))

       do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         wallThickness(fl) = thick_wall(l)
         roofThickness(fl) = thick_roof(l)
       end do

       call UrbanSetBuildingWallThickness(urban, c_loc(wallThickness), &
            num_urbanl, status)
       if (status /= URBAN_SUCCESS) then
         write(iulog,*) 'ERROR: UrbanSetBuildingWallThickness failed with status: ', status
         call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       call UrbanSetBuildingRoofThickness(urban, c_loc(roofThickness), &
            num_urbanl, status)
       if (status /= URBAN_SUCCESS) then
         write(iulog,*) 'ERROR: UrbanSetBuildingRoofThickness failed with status: ', status
         call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       deallocate(wallThickness)
       deallocate(roofThickness)

     end associate

   end subroutine SetBuildingTemperature

   !-----------------------------------------------------------------------
   subroutine SetUrbanParameters(urban, num_urbanl, filter_urbanl, &
        urbanparams_vars, frictionvel_vars, soilstate_vars)
     !
     implicit none
     !
     type(UrbanType)        , intent(inout) :: urban
     integer(c_int)         , intent(in)    :: num_urbanl
     integer                , intent(in)    :: filter_urbanl(:) ! urban landunit filter
     type(urbanparams_type) , intent(in)    :: urbanparams_vars
     type(frictionvel_type) , intent(in)    :: frictionvel_vars
     type(soilstate_type)   , intent(in)    :: soilstate_vars

     call SetCanyonHwr(urban, num_urbanl, filter_urbanl)
     call SetFracPervRoadOfTotalRoad(urban, num_urbanl, filter_urbanl)
     call SetWtRoof(urban, num_urbanl, filter_urbanl)
     call SetHeightParameters(urban, num_urbanl, filter_urbanl, &
          urbanparams_vars, frictionvel_vars)
     call SetAlbedo(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetEmissivity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetNumberOfActiveLayersImperviousRoad(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetBuildingTemperature(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetThermalConductivity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetHeatCapacity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetSoilProperties(urban, num_urbanl, filter_urbanl, soilstate_vars)
     call SetCanyonAirStates(urban, num_urbanl, filter_urbanl)

   end subroutine SetUrbanParameters

end module UrbanxxMod
