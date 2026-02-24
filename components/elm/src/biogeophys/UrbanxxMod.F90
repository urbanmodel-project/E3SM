module UrbanxxMod

  use iso_c_binding
  use urban_mod
  use urban_kokkos_interface
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use decompMod            , only : bounds_type
  use UrbanParamsType      , only : urbanparams_type
  use SolarAbsorbedType    , only : solarabs_type
  use SurfaceAlbedoType    , only : surfalb_type
  use LandunitType         , only : lun_pp
  use FrictionVelocityType , only : frictionvel_type
  use ColumnDataType       , only : col_es, col_pp
  use LandunitDataType     , only : lun_es, lun_ws
  use SoilStateType        , only : soilstate_type
  use abortutils           , only : endrun
  use UrbanxxInstanceMod   , only : urbanxx, numBands, numTypes

  implicit none

  private

  public :: urbanxx_initialize
  public :: SetHeightParameters
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

    call SetUrbanParameters(urbanxx, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, &
         urbanparams_vars, frictionvel_vars, soilstate_vars)

    ! Setup urban model (initialize temperatures and other setup tasks)
    call UrbanSetup(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    if (masterproc) then
       write(*,*) 'Completed urban model setup'
    end if

    ! Allocate persistent buffers for physics modules
    ! (init calls moved to elm_initializeMod.F90 to avoid circular dependency)

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
      do j = 1, nlevurb
         do fl = 1, num_urbanl
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
   subroutine SetSoilProperties(urban, num_urbanl, num_urbanc, filter_urbanc, soilstate_vars)
     !
     use elm_varpar, only : nlevgrnd
     use ColumnType, only : col_pp
     use column_varcon, only : icol_road_perv
     !
     implicit none
     !
     type(UrbanType)      , intent(in)    :: urban
     integer(c_int)       , intent(in)    :: num_urbanl
     integer(c_int)       , intent(in)    :: num_urbanc
     integer              , intent(in)    :: filter_urbanc(:) ! urban column filter
     type(soilstate_type) , intent(in)    :: soilstate_vars
     !
     integer(c_int)                       :: status
     integer                              :: fc, c, j, idx, nlevbed
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
          nlev2bed => col_pp%nlevbed                & ! Input: [integer (:)] number of layers to bedrock
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
           do fc = 1, num_urbanc
             c = filter_urbanc(fc)
             if (col_pp%itype(c) == icol_road_perv) then
               idx = idx + 1
               nlevbed = nlev2bed(c)
               if (j <= nlevbed) then
                 sand(idx) = cellsand(c, j)
                 clay(idx) = cellclay(c, j)
                 organic(idx) = cellorg(c, j)
               else
                 ! Below bedrock: set sand to 100%, clay to 0%
                 sand(idx) = 100.0_r8
                 clay(idx) = 0.0_r8
                 organic(idx) = 0.0_r8
               end if
             end if
           end do
         end do
       else
         ! LayoutRight: Last dimension (layers) varies fastest
         ! Iterate: landunits (outer), layer (inner)
         idx = 0
         do fc = 1, num_urbanc
           c = filter_urbanc(fc)
           if (col_pp%itype(c) == icol_road_perv) then
             nlevbed = nlev2bed(c)
             do j = 1, nlevgrnd
               idx = idx + 1
               if (j <= nlevbed) then
                 sand(idx) = cellsand(c, j)
                 clay(idx) = cellclay(c, j)
                 organic(idx) = cellorg(c, j)
               else
                 ! Below bedrock: set sand to 100%, clay to 0%
                 sand(idx) = 100.0_r8
                 clay(idx) = 0.0_r8
                 organic(idx) = 0.0_r8
               end if
             end do
           end if
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
   subroutine SetLayerTemperatures(urban, num_urbanl, filter_urbanl)
     !
     use elm_varpar, only : nlevurb, nlevgrnd
     use column_varcon, only : icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, icol_shadewall
     !
     implicit none
     !
     type(UrbanType)        , intent(in) :: urban
     integer(c_int)         , intent(in) :: num_urbanl
     integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
     !
     integer(c_int)                       :: status
     integer                              :: fl, l, c, c_start, c_end, j
     integer                              :: idx_roof, idx_improad, idx_pervroad, idx_sunwall, idx_shadwall
     integer(c_int), dimension(2)         :: size2D_urban, size2D_soil
     logical(c_bool)                      :: isLayoutLeft
     real(c_double), allocatable, target  :: tempRoof(:)
     real(c_double), allocatable, target  :: tempImpRoad(:)
     real(c_double), allocatable, target  :: tempPervRoad(:)
     real(c_double), allocatable, target  :: tempSunlitWall(:)
     real(c_double), allocatable, target  :: tempShadedWall(:)

     associate(                       &
          ctype  =>    col_pp%itype , & ! Input:  [integer (:)    ]  column type
          coli   =>    lun_pp%coli  , & ! Input:  [integer (:)    ]  beginning column index for landunit
          colf   =>    lun_pp%colf  , & ! Input:  [integer (:)    ]  ending column index for landunit
          t_soisno =>  col_es%t_soisno & ! Input:  [real(r8) (:,:) ]  soil temperature (K)
          )

       ! Check memory layout
       isLayoutLeft = UrbanKokkosIsLayoutLeft()

       ! Allocate arrays for layer temperatures
       ! Urban surfaces: (num_urbanl, nlevurb)
       allocate(tempRoof(num_urbanl * nlevurb))
       allocate(tempSunlitWall(num_urbanl * nlevurb))
       allocate(tempShadedWall(num_urbanl * nlevurb))
       
       ! Road surfaces: (num_urbanl, nlevgrnd)
       allocate(tempImpRoad(num_urbanl * nlevgrnd))
       allocate(tempPervRoad(num_urbanl * nlevgrnd))

       ! Extract layer temperatures from columns based on column type
       if (isLayoutLeft) then
         ! LayoutLeft: iterate layers in outer loop, landunits in inner loop
         ! For urban surfaces (roof, walls)
         idx_roof = 0
         idx_sunwall = 0
         idx_shadwall = 0
         do j = 1, nlevurb
           do fl = 1, num_urbanl
             l = filter_urbanl(fl)
             c_start = coli(l)
             c_end   = colf(l)
             
             do c = c_start, c_end
               select case (ctype(c))
               case (icol_roof)
                 idx_roof = idx_roof + 1
                 tempRoof(idx_roof) = t_soisno(c, j)
               case (icol_sunwall)
                 idx_sunwall = idx_sunwall + 1
                 tempSunlitWall(idx_sunwall) = t_soisno(c, j)
               case (icol_shadewall)
                 idx_shadwall = idx_shadwall + 1
                 tempShadedWall(idx_shadwall) = t_soisno(c, j)
               end select
             end do
           end do
         end do
         
         ! For road surfaces
         idx_improad = 0
         idx_pervroad = 0
         do j = 1, nlevgrnd
           do fl = 1, num_urbanl
             l = filter_urbanl(fl)
             c_start = coli(l)
             c_end   = colf(l)
             
             do c = c_start, c_end
               select case (ctype(c))
               case (icol_road_imperv)
                 idx_improad = idx_improad + 1
                 tempImpRoad(idx_improad) = t_soisno(c, j)
               case (icol_road_perv)
                 idx_pervroad = idx_pervroad + 1
                 tempPervRoad(idx_pervroad) = t_soisno(c, j)
               end select
             end do
           end do
         end do
       else
         ! LayoutRight: iterate landunits in outer loop, layers in inner loop
         ! For urban surfaces (roof, walls)
         idx_roof = 0
         idx_sunwall = 0
         idx_shadwall = 0
         do fl = 1, num_urbanl
           l = filter_urbanl(fl)
           c_start = coli(l)
           c_end   = colf(l)
           
           do c = c_start, c_end
             if (ctype(c) == icol_roof) then
               do j = 1, nlevurb
                 idx_roof = idx_roof + 1
                 tempRoof(idx_roof) = t_soisno(c, j)
               end do
             else if (ctype(c) == icol_sunwall) then
               do j = 1, nlevurb
                 idx_sunwall = idx_sunwall + 1
                 tempSunlitWall(idx_sunwall) = t_soisno(c, j)
               end do
             else if (ctype(c) == icol_shadewall) then
               do j = 1, nlevurb
                 idx_shadwall = idx_shadwall + 1
                 tempShadedWall(idx_shadwall) = t_soisno(c, j)
               end do
             end if
           end do
         end do
         
         ! For road surfaces
         idx_improad = 0
         idx_pervroad = 0
         do fl = 1, num_urbanl
           l = filter_urbanl(fl)
           c_start = coli(l)
           c_end   = colf(l)
           
           do c = c_start, c_end
             if (ctype(c) == icol_road_imperv) then
               do j = 1, nlevgrnd
                 idx_improad = idx_improad + 1
                 tempImpRoad(idx_improad) = t_soisno(c, j)
               end do
             else if (ctype(c) == icol_road_perv) then
               do j = 1, nlevgrnd
                 idx_pervroad = idx_pervroad + 1
                 tempPervRoad(idx_pervroad) = t_soisno(c, j)
               end do
             end if
           end do
         end do
       end if

       ! Set size arrays
       size2D_urban(1) = num_urbanl
       size2D_urban(2) = nlevurb
       size2D_soil(1) = num_urbanl
       size2D_soil(2) = nlevgrnd

       ! Set roof layer temperatures
       call UrbanSetLayerTempRoof(urban, c_loc(tempRoof), size2D_urban, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       ! Set impervious road layer temperatures
       call UrbanSetLayerTempImperviousRoad(urban, c_loc(tempImpRoad), &
         size2D_soil, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       ! Set pervious road layer temperatures
       call UrbanSetLayerTempPerviousRoad(urban, c_loc(tempPervRoad), &
         size2D_soil, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       ! Set sunlit wall layer temperatures
       call UrbanSetLayerTempSunlitWall(urban, c_loc(tempSunlitWall), &
         size2D_urban, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       ! Set shaded wall layer temperatures
       call UrbanSetLayerTempShadedWall(urban, c_loc(tempShadedWall), &
         size2D_urban, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       if (masterproc) then
         write(iulog,*) 'Set layer temperatures for all urban surfaces from ELM data'
       end if

       deallocate(tempRoof)
       deallocate(tempImpRoad)
       deallocate(tempPervRoad)
       deallocate(tempSunlitWall)
       deallocate(tempShadedWall)
     end associate

   end subroutine SetLayerTemperatures

   !-----------------------------------------------------------------------
   subroutine SetUrbanParameters(urban, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, &
        urbanparams_vars, frictionvel_vars, soilstate_vars)
     !
     implicit none
     !
     type(UrbanType)        , intent(inout) :: urban
     integer(c_int)         , intent(in)    :: num_urbanl
     integer                , intent(in)    :: filter_urbanl(:) ! urban landunit filter
     integer(c_int)         , intent(in)    :: num_urbanc
     integer                , intent(in)    :: filter_urbanc(:) ! urban column filter
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
     call SetLayerTemperatures(urban, num_urbanl, filter_urbanl)
     call SetThermalConductivity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetHeatCapacity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetSoilProperties(urban, num_urbanl, num_urbanc, filter_urbanc, soilstate_vars)
     call SetCanyonAirStates(urban, num_urbanl, filter_urbanl)

   end subroutine SetUrbanParameters

end module UrbanxxMod
