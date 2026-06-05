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
  use SoilHydrologyType    , only : soilhydrology_type
  use abortutils           , only : endrun
  use UrbanxxInstanceMod   , only : urbanxx, numBands, numTypes
  use TopounitDataType     , only : topounit_atmospheric_state

  implicit none

  ! Set to .true. to override URBANxx-internal pedotransfer results with ELM values.
  ! Avoids precision differences from Kokkos::pow vs Fortran ** in the same formulas.
  logical, parameter :: use_elm_pervroad_derived_props = .true.

  private

  public :: urbanxx_initialize
  public :: SetHeightParameters
contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_initialize(bounds, num_urbanl, filter_urbanl, &
       num_urbanc, filter_urbanc, num_urbanp, filter_urbanp, &
       urbanparams_vars, solarabs_vars, surfalb_vars, top_as, &
       soilstate_vars, soilhydrology_vars)
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
    type(topounit_atmospheric_state) , intent(in) :: top_as
    type(soilstate_type)   , intent(in) :: soilstate_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars

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
         urbanparams_vars, top_as, soilstate_vars, soilhydrology_vars)

    ! Setup urban model (initialize temperatures and other setup tasks)
    call UrbanSetup(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
    if (masterproc) then
       write(*,*) 'Completed urban model setup'
    end if

    ! UrbanSetup (above) calls UrbanInitializePerviousRoadSoils which resets
    ! perviousRoad WaterLiquid/WaterIce to a default value.  Re-seed them from
    ! ELM's h2osoi_liq/ice so that thermal conductivity and heat capacity use
    ! the correct moisture on the first post-restart timestep.
    call SetPerviousRoadSoilWater(urbanxx, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)
    if (masterproc) then
       write(iulog,*) 'Re-seeded pervious road WaterLiquid/Ice after UrbanSetup from ELM data'
    end if

    ! Seed water table depth (Zwt) and unconfined aquifer water (Wa) from ELM
    ! state so that URBANxx hydrology starts from the correct state on restart.
    ! UrbanInitializePerviousRoadSoils sets Zwt=4.8 m and Wa=0 by default.
    call SetWaterTableFromState(urbanxx, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, &
         soilhydrology_vars)
    if (masterproc) then
       write(iulog,*) 'Re-seeded pervious road Zwt/Wa from ELM data'
    end if

    ! Override derived soil properties with ELM values to eliminate pedotransfer
    ! precision differences (Kokkos::pow vs Fortran **).
    if (use_elm_pervroad_derived_props) then
       if (masterproc) then
          write(*,*) 'Set pervious road soil properties derived from ELM values'
       end if
       call SetPerviousRoadDerivedProperties(urbanxx, num_urbanl, num_urbanc, &
            filter_urbanc, soilstate_vars)
    end if

    ! Seed building temperature from the restored ELM restart state so that
    ! the HVAC boundary condition is correct on the first post-restart timestep.
    call SetBuildingTemperatureFromState(urbanxx, num_urbanl, filter_urbanl)

    ! Seed AbsorbedShortRad (absRoof, absImpRoad, absPerRoad, absSunWall,
    ! absShadWall) from the ELM restart state.  These views are used by
    ! UrbanComputeNetShortwaveRadiation (called every timestep via
    ! urbanxx_SetAtmosphericForcing) before UrbanComputeNetShortwave (step 5)
    ! has had a chance to recompute them.  Without seeding they are zero on the
    ! first post-restart timestep, causing swnet_roof = 0 and a large top-BC
    ! error in the roof heat diffusion.
    call SetShortwaveAlbedo(urbanxx, num_urbanl, filter_urbanl, solarabs_vars)
    if (masterproc) then
       write(iulog,*) 'Re-seeded shortwave albedo (AbsorbedShortRad) from ELM solarabs_vars restart data'
    end if

    ! Seed EFluxForAC (previous-timestep AC heat) from ELM's eflx_urban_ac so
    ! that the road boundary condition in UrbanComputeHeatDiffusion is correct
    ! on the first post-restart timestep.  On cold start eflx_urban_ac = 0.
    call SetACHeatFromState(urbanxx, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)

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
       urbanparams_vars, top_as)
    !
    implicit none
    !
    type(UrbanType)        , intent(in) :: urban
    integer(c_int)         , intent(in) :: num_urbanl
    integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    type(topounit_atmospheric_state) , intent(in) :: top_as
    !
    integer(c_int)                       :: status
    integer                              :: fl, fp, l, p, t
    real(c_double) , allocatable, target :: forcHgtT(:)
    real(c_double) , allocatable, target :: forcHgtU(:)
    real(c_double) , allocatable, target :: zDTown(:)
    real(c_double) , allocatable, target :: z0Town(:)
    real(c_double) , allocatable, target :: htRoof(:)
    real(c_double) , allocatable, target :: windHgtCanyon(:)

    associate(                                                          &
         forc_hgt_t       => top_as%zbot                              , & ! Input:  [real(r8) (:)   ] observational height of temperature [m]
         forc_hgt_u       => top_as%zbot                              , & ! Input:  [real(r8) (:)   ] observational height of wind [m]
         forc_hgt_q       => top_as%zbot                              , & ! Input:  [real(r8) (:)   ] observational height of specific humidity [m]
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
         t = lun_pp%topounit(l)
         forcHgtT(fl)      = forc_hgt_t(t) + z_0_town(l) + z_d_town(l)
         forcHgtU(fl)      = forc_hgt_u(t) + z_0_town(l) + z_d_town(l)
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

      !if (masterproc) then
      !   write(*,*) 'Set height parameters'
      !end if

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
                 ! Below bedrock: ELM duplicates last active soil layer sand/clay to
                 ! compute watsat for bedrock layers while overriding csol=csol_bedrock.
                 ! Match that here so URBANxx watsat at bedrock equals ELM's watsat.
                 sand(idx) = cellsand(c, nlevbed)
                 clay(idx) = cellclay(c, nlevbed)
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
                 ! Below bedrock: duplicate last active layer texture (see ELM SoilStateType.F90)
                 sand(idx) = cellsand(c, nlevbed)
                 clay(idx) = cellclay(c, nlevbed)
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
   subroutine SetPerviousRoadDerivedProperties(urban, num_urbanl, num_urbanc, &
       filter_urbanc, soilstate_vars)
     !
     use elm_varpar,    only : nlevgrnd
     use ColumnType,    only : col_pp
     use column_varcon, only : icol_road_perv
     !
     implicit none
     !
     type(UrbanType)      , intent(in) :: urban
     integer(c_int)       , intent(in) :: num_urbanl
     integer(c_int)       , intent(in) :: num_urbanc
     integer              , intent(in) :: filter_urbanc(:)
     type(soilstate_type) , intent(in) :: soilstate_vars
     !
     integer(c_int)                       :: status
     integer                              :: fc, c, j, idx, nlevbed
     integer(c_int)                       :: totalSize
     integer(c_int), dimension(2)         :: size2D
     logical(c_bool)                      :: isLayoutLeft
     real(c_double), allocatable, target  :: watsat(:)
     real(c_double), allocatable, target  :: watdry(:)
     real(c_double), allocatable, target  :: watopt(:)
     real(c_double), allocatable, target  :: bsw(:)
     real(c_double), allocatable, target  :: sucsat(:)
     real(c_double), allocatable, target  :: hksat(:)
     real(c_double), allocatable, target  :: tkdry(:)
     real(c_double), allocatable, target  :: tksat(:)
     real(c_double), allocatable, target  :: tkminerals(:)
     real(c_double), allocatable, target  :: cvsolids(:)

     associate(                                   &
          watsat_col => soilstate_vars%watsat_col, & ! Input: [real(r8) (:,:)] volumetric soil water at saturation
          watdry_col => soilstate_vars%watdry_col, & ! Input: [real(r8) (:,:)] vol. water at wilting point (btran=0)
          watopt_col => soilstate_vars%watopt_col, & ! Input: [real(r8) (:,:)] vol. water at field-capacity analogue (btran=1)
          bsw_col    => soilstate_vars%bsw_col   , & ! Input: [real(r8) (:,:)] Clapp-Hornberger parameter b
          sucsat_col => soilstate_vars%sucsat_col , & ! Input: [real(r8) (:,:)] minimum soil suction [mm]
          hksat_col  => soilstate_vars%hksat_col  , & ! Input: [real(r8) (:,:)] hydraulic conductivity at saturation [mm/s]
          tkdry_col  => soilstate_vars%tkdry_col  , & ! Input: [real(r8) (:,:)] thermal conductivity of dry soil [W/(m K)]
          tksatu_col => soilstate_vars%tksatu_col , & ! Input: [real(r8) (:,:)] thermal conductivity, saturated soil [W/(m K)]
          tkmg_col   => soilstate_vars%tkmg_col   , & ! Input: [real(r8) (:,:)] thermal conductivity of soil minerals [W/(m K)]
          csol_col   => soilstate_vars%csol_col   , & ! Input: [real(r8) (:,:)] heat capacity of soil solids [J/(m3 K)]
          nlev2bed   => col_pp%nlevbed              & ! Input: [integer (:)] number of layers to bedrock
          )

       totalSize = num_urbanl * nlevgrnd
       size2D(1) = num_urbanl
       size2D(2) = nlevgrnd

       allocate(watsat(totalSize))
       allocate(watdry(totalSize))
       allocate(watopt(totalSize))
       allocate(bsw(totalSize))
       allocate(sucsat(totalSize))
       allocate(hksat(totalSize))
       allocate(tkdry(totalSize))
       allocate(tksat(totalSize))
       allocate(tkminerals(totalSize))
       allocate(cvsolids(totalSize))

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
               tkdry(idx)     = tkdry_col(c, j)
               tksat(idx)     = tksatu_col(c, j)
               tkminerals(idx)= tkmg_col(c, j)
               cvsolids(idx)  = csol_col(c, j)
               if (j <= nlevbed) then
                 watsat(idx)    = watsat_col(c, j)
                 watdry(idx)    = watdry_col(c, j)
                 watopt(idx)    = watopt_col(c, j)
                 bsw(idx)       = bsw_col(c, j)
                 sucsat(idx)    = sucsat_col(c, j)
                 hksat(idx)     = hksat_col(c, j)
               else
                 ! Below bedrock: duplicate last active layer
                 watsat(idx)    = watsat_col(c, nlevbed)
                 watdry(idx)    = watdry_col(c, nlevbed)
                 watopt(idx)    = watopt_col(c, nlevbed)
                 bsw(idx)       = bsw_col(c, nlevbed)
                 sucsat(idx)    = sucsat_col(c, nlevbed)
                 hksat(idx)     = hksat_col(c, nlevbed)
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
               tkdry(idx)     = tkdry_col(c, j)
               tksat(idx)     = tksatu_col(c, j)
               tkminerals(idx)= tkmg_col(c, j)
               cvsolids(idx)  = csol_col(c, j)
               if (j <= nlevbed) then
                 watsat(idx)    = watsat_col(c, j)
                 watdry(idx)    = watdry_col(c, j)
                 watopt(idx)    = watopt_col(c, j)
                 bsw(idx)       = bsw_col(c, j)
                 sucsat(idx)    = sucsat_col(c, j)
                 hksat(idx)     = hksat_col(c, j)
               else
                 ! Below bedrock: duplicate last active layer
                 watsat(idx)    = watsat_col(c, nlevbed)
                 watdry(idx)    = watdry_col(c, nlevbed)
                 watopt(idx)    = watopt_col(c, nlevbed)
                 bsw(idx)       = bsw_col(c, nlevbed)
                 sucsat(idx)    = sucsat_col(c, nlevbed)
                 hksat(idx)     = hksat_col(c, nlevbed)
               end if
             end do
           end if
         end do
       end if

       call UrbanSetWatSatForPerviousRoad(urban, c_loc(watsat), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetWatDryForPerviousRoad(urban, c_loc(watdry), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetWatOptForPerviousRoad(urban, c_loc(watopt), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetBswForPerviousRoad(urban, c_loc(bsw), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetSucSatForPerviousRoad(urban, c_loc(sucsat), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetHkSatForPerviousRoad(urban, c_loc(hksat), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetTkDryForPerviousRoad(urban, c_loc(tkdry), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetTkSatForPerviousRoad(urban, c_loc(tksat), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetTkMineralsForPerviousRoad(urban, c_loc(tkminerals), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetCvSolidsForPerviousRoad(urban, c_loc(cvsolids), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       if (masterproc) then
         write(iulog,*) 'Set derived soil properties for pervious road from ELM data (watsat,watdry,watopt,bsw,sucsat,hksat,tkdry,tksat,tkminerals,cvsolids)'
       end if

       deallocate(watsat)
       deallocate(watdry)
       deallocate(watopt)
       deallocate(bsw)
       deallocate(sucsat)
       deallocate(hksat)
       deallocate(tkdry)
       deallocate(tksat)
       deallocate(tkminerals)
       deallocate(cvsolids)
     end associate

   end subroutine SetPerviousRoadDerivedProperties

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
   subroutine SetBuildingTemperatureFromState(urban, num_urbanl, filter_urbanl)
     !
     ! !DESCRIPTION:
     ! Seed URBANxx BuildingTemperature from ELM's t_building_lun.  Called once
     ! during urbanxx_initialize so that restart runs pick up the saved building
     ! temperature rather than the cold-start value from UrbanSetup.
     ! Uses module-level lun_es (same pattern as col_es, col_ws, etc.)
     ! On cold start lun_es%t_building is spval — skip seeding in that case to
     ! preserve the physically initialised value from UrbanSetup.
     !
     use elm_varcon , only : spval
     implicit none
     !
     type(UrbanType)        , intent(in) :: urban
     integer(c_int)         , intent(in) :: num_urbanl
     integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
     !
     integer(c_int)                       :: status
     integer                              :: fl, l
     real(c_double) , allocatable, target :: tBuilding(:)

     allocate(tBuilding(num_urbanl))

     ! On a cold start t_building is spval (not yet computed by UrbanFluxes).
     ! Only seed URBANxx when ELM has a real value (e.g. from a restart file).
     if (lun_es%t_building(filter_urbanl(1)) /= spval) then
       do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         tBuilding(fl) = real(lun_es%t_building(l), c_double)
       end do

       call UrbanSetBuildingTemperature(urban, c_loc(tBuilding), num_urbanl, status)
       if (status /= URBAN_SUCCESS) then
         write(iulog,*) 'ERROR: UrbanSetBuildingTemperature failed with status: ', status
         call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
     end if

     deallocate(tBuilding)

   end subroutine SetBuildingTemperatureFromState

   !-----------------------------------------------------------------------
   subroutine SetShortwaveAlbedo(urban, num_urbanl, filter_urbanl, solarabs_vars)
     !
     ! !DESCRIPTION:
     ! Seed URBANxx's AbsorbedShortRad views (absRoof, absImpRoad, absPerRoad,
     ! absSunWall, absShadWall) from ELM's sabs_*_dir_lun / sabs_*_dif_lun.
     ! Despite the ELM variable names, these quantities are shortwave albedos
     ! (absorbed radiation per unit incoming radiation), not absolute fluxes.
     ! Called once during urbanxx_initialize so that restart runs have non-zero
     ! absorbed shortwave fractions before UrbanComputeNetShortwave (step 5)
     ! has run.  Without seeding, UrbanComputeNetShortwaveRadiation (called
     ! every timestep) computes NetShortRad = abs * ForcSRad = 0 on the first
     ! post-restart timestep, producing a large top-BC error in heat diffusion.
     ! Indexing: itype=0 -> DIFFUSE (*_dif), itype=1 -> DIRECT (*_dir)
     ! On cold start sabs_*_lun = 0, so seeding is always safe.
     !
     use elm_varcon , only : spval
     implicit none
     !
     type(UrbanType)        , intent(in) :: urban
     integer(c_int)         , intent(in) :: num_urbanl
     integer                , intent(in) :: filter_urbanl(:) ! urban landunit filter
     type(solarabs_type)    , intent(in) :: solarabs_vars
     !
     integer(c_int)                       :: status
     integer                              :: fl, l, iband, itype, count
     integer(c_int)                       :: totalSize3D
     integer(c_int), dimension(3)         :: size3D
     real(c_double) , allocatable, target :: absRoof(:)
     real(c_double) , allocatable, target :: absImpRoad(:)
     real(c_double) , allocatable, target :: absPerRoad(:)
     real(c_double) , allocatable, target :: absSunWall(:)
     real(c_double) , allocatable, target :: absShadWall(:)

     associate(                                                              &
          sabs_roof_dir      => solarabs_vars%sabs_roof_dir_lun          , & ! Input: [real(r8) (:,:)]
          sabs_roof_dif      => solarabs_vars%sabs_roof_dif_lun          , & ! Input: [real(r8) (:,:)]
          sabs_improad_dir   => solarabs_vars%sabs_improad_dir_lun       , & ! Input: [real(r8) (:,:)]
          sabs_improad_dif   => solarabs_vars%sabs_improad_dif_lun       , & ! Input: [real(r8) (:,:)]
          sabs_perroad_dir   => solarabs_vars%sabs_perroad_dir_lun       , & ! Input: [real(r8) (:,:)]
          sabs_perroad_dif   => solarabs_vars%sabs_perroad_dif_lun       , & ! Input: [real(r8) (:,:)]
          sabs_sunwall_dir   => solarabs_vars%sabs_sunwall_dir_lun       , & ! Input: [real(r8) (:,:)]
          sabs_sunwall_dif   => solarabs_vars%sabs_sunwall_dif_lun       , & ! Input: [real(r8) (:,:)]
          sabs_shadewall_dir => solarabs_vars%sabs_shadewall_dir_lun     , & ! Input: [real(r8) (:,:)]
          sabs_shadewall_dif => solarabs_vars%sabs_shadewall_dif_lun       & ! Input: [real(r8) (:,:)]
          )

       size3D = [num_urbanl, numBands, numTypes]
       totalSize3D = num_urbanl * numBands * numTypes

       allocate(absRoof(totalSize3D))
       allocate(absImpRoad(totalSize3D))
       allocate(absPerRoad(totalSize3D))
       allocate(absSunWall(totalSize3D))
       allocate(absShadWall(totalSize3D))

       if (UrbanKokkosIsLayoutLeft()) then

          ! Kokkos layout left: 'l'-index incremented fastest
          count = 0
          do itype = 0, 1
             do iband = 0, numBands - 1
                do fl = 1, num_urbanl
                   l = filter_urbanl(fl)
                   count = count + 1
                   if (itype == 1) then
                      ! itype = 1: diffuse
                      absRoof(count)    = real(sabs_roof_dif(l, iband+1),      c_double)
                      absImpRoad(count) = real(sabs_improad_dif(l, iband+1),   c_double)
                      absPerRoad(count) = real(sabs_perroad_dif(l, iband+1),   c_double)
                      absSunWall(count) = real(sabs_sunwall_dif(l, iband+1),   c_double)
                      absShadWall(count)= real(sabs_shadewall_dif(l, iband+1), c_double)
                   else
                      ! itype = 0: direct
                      absRoof(count)    = real(sabs_roof_dir(l, iband+1),      c_double)
                      absImpRoad(count) = real(sabs_improad_dir(l, iband+1),   c_double)
                      absPerRoad(count) = real(sabs_perroad_dir(l, iband+1),   c_double)
                      absSunWall(count) = real(sabs_sunwall_dir(l, iband+1),   c_double)
                      absShadWall(count)= real(sabs_shadewall_dir(l, iband+1), c_double)
                   end if
                end do
             end do
          end do

       else

          ! Kokkos layout right: 'itype'-index incremented fastest
          count = 0
          do fl = 1, num_urbanl
             l = filter_urbanl(fl)
             do iband = 0, numBands - 1
                do itype = 0, 1
                   count = count + 1
                   if (itype == 1) then
                      ! itype = 1: diffuse
                      absRoof(count)    = real(sabs_roof_dif(l, iband+1),      c_double)
                      absImpRoad(count) = real(sabs_improad_dif(l, iband+1),   c_double)
                      absPerRoad(count) = real(sabs_perroad_dif(l, iband+1),   c_double)
                      absSunWall(count) = real(sabs_sunwall_dif(l, iband+1),   c_double)
                      absShadWall(count)= real(sabs_shadewall_dif(l, iband+1), c_double)
                   else
                      ! itype = 0: direct
                      absRoof(count)    = real(sabs_roof_dir(l, iband+1),      c_double)
                      absImpRoad(count) = real(sabs_improad_dir(l, iband+1),   c_double)
                      absPerRoad(count) = real(sabs_perroad_dir(l, iband+1),   c_double)
                      absSunWall(count) = real(sabs_sunwall_dir(l, iband+1),   c_double)
                      absShadWall(count)= real(sabs_shadewall_dir(l, iband+1), c_double)
                   end if
                end do
             end do
          end do

       end if

       call UrbanSetAbsorbedShortRadRoof(urban, c_loc(absRoof), size3D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetAbsorbedShortRadImperviousRoad(urban, c_loc(absImpRoad), size3D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetAbsorbedShortRadPerviousRoad(urban, c_loc(absPerRoad), size3D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetAbsorbedShortRadSunlitWall(urban, c_loc(absSunWall), size3D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetAbsorbedShortRadShadedWall(urban, c_loc(absShadWall), size3D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       deallocate(absRoof)
       deallocate(absImpRoad)
       deallocate(absPerRoad)
       deallocate(absSunWall)
       deallocate(absShadWall)

     end associate

   end subroutine SetShortwaveAlbedo

   !-----------------------------------------------------------------------
   subroutine SetACHeatFromState(urban, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)
     !
     ! !DESCRIPTION:
     ! Seed URBANxx's EFluxForAC from ELM's eflx_urban_ac for each urban
     ! landunit.  Computes eflx_heat_from_ac_patch (the per-road-area AC heat
     ! flux) using the same weighted formula as ELM's UrbanFluxesMod:
     !   eflx_heat_from_ac(l) = wt_roof*|ac_roof| + (1-wt_roof)*hwr*(|ac_sunwall|+|ac_shadwall|)
     !   eflx_heat_from_ac_patch = eflx_heat_from_ac(l) / (1 - wt_roof)
     ! Called once during urbanxx_initialize so that restart runs start with
     ! the correct previous-timestep AC heat used as the road boundary condition
     ! in UrbanComputeHeatDiffusion.
     ! On cold start eflx_urban_ac(c) = 0, so seeding is always safe.
     !
     use ColumnType     , only : col_pp
     use ColumnDataType , only : col_ef
     use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
     use LandunitType   , only : lun_pp
     !
     implicit none
     !
     type(UrbanType)        , intent(in) :: urban
     integer(c_int)         , intent(in) :: num_urbanl
     integer                , intent(in) :: filter_urbanl(:)
     integer(c_int)         , intent(in) :: num_urbanc
     integer                , intent(in) :: filter_urbanc(:)
     !
     integer(c_int)                       :: status
     integer                              :: fl, l, fc, c
     real(c_double), allocatable, target  :: acFlux(:)
     real(r8) :: ac_roof, ac_sunwall, ac_shadwall, wt_roof, hwr, eflx_heat_from_ac_l

     associate( &
         eflx_urban_ac => col_ef%eflx_urban_ac  , &
         canyon_hwr    => lun_pp%canyon_hwr      , &
         wtlunit_roof  => lun_pp%wtlunit_roof      &
         )

       allocate(acFlux(num_urbanl))
       acFlux = 0._c_double

       ! For each urban landunit, compute eflx_heat_from_ac_patch:
       !   eflx_heat_from_ac_patch = [wt_roof*|ac_roof| + (1-wt_roof)*hwr*(|ac_sunwall|+|ac_shadwall|)]
       !                             / (1 - wt_roof)
       do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         wt_roof    = wtlunit_roof(l)
         hwr        = canyon_hwr(l)
         ac_roof    = 0._r8
         ac_sunwall = 0._r8
         ac_shadwall= 0._r8
         do fc = 1, num_urbanc
           c = filter_urbanc(fc)
           if (col_pp%landunit(c) /= l) cycle
           select case (col_pp%itype(c))
           case (icol_roof)
             ac_roof     = abs(eflx_urban_ac(c))
           case (icol_sunwall)
             ac_sunwall  = abs(eflx_urban_ac(c))
           case (icol_shadewall)
             ac_shadwall = abs(eflx_urban_ac(c))
           end select
         end do
         eflx_heat_from_ac_l = wt_roof * ac_roof + (1._r8 - wt_roof) * hwr * (ac_sunwall + ac_shadwall)
         if ((1._r8 - wt_roof) > 0._r8) then
           acFlux(fl) = real(eflx_heat_from_ac_l / (1._r8 - wt_roof), c_double)
         end if
       end do

       call UrbanSetEFluxForAC(urban, c_loc(acFlux), num_urbanl, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       deallocate(acFlux)

     end associate

     if (masterproc) then
       write(iulog,*) 'Seeded EFluxForAC (AC heat) for all urban landunits from ELM data'
     end if

   end subroutine SetACHeatFromState

   !-----------------------------------------------------------------------
   subroutine SetImperviousRoadSoilWater(urban, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)
     !
     ! !DESCRIPTION:
     ! Initialize URBANxx's persistent WaterLiquid/WaterIce for impervious road
     ! deep soil layers from ELM's initial h2osoi_liq/h2osoi_ice.  Called once
     ! at startup/restart alongside SetLayerTemperatures.  After initialization,
     ! ComputeHeatDiffusion maintains these arrays internally via phase change.
     !
     use elm_varpar     , only : nlevgrnd
     use ColumnType     , only : col_pp
     use ColumnDataType , only : col_ws
     use column_varcon  , only : icol_road_imperv
     use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
     !
     implicit none
     !
     type(UrbanType)  , intent(in) :: urban
     integer(c_int)   , intent(in) :: num_urbanl
     integer          , intent(in) :: filter_urbanl(:)
     integer(c_int)   , intent(in) :: num_urbanc
     integer          , intent(in) :: filter_urbanc(:)
     !
     integer(c_int)                       :: status
     integer                              :: fc, c, j, idx
     integer(c_int), dimension(2)         :: size2D
     logical(c_bool)                      :: isLayoutLeft
     real(c_double), allocatable, target  :: h2oLiq(:)
     real(c_double), allocatable, target  :: h2oIce(:)

     associate( &
         h2osoi_liq => col_ws%h2osoi_liq, & ! Input: [real(r8) (:,:) ] liquid water (kg/m2)
         h2osoi_ice => col_ws%h2osoi_ice  & ! Input: [real(r8) (:,:) ] ice lens (kg/m2)
         )

       allocate(h2oLiq(num_urbanl * nlevgrnd))
       allocate(h2oIce(num_urbanl * nlevgrnd))
       size2D(1) = num_urbanl
       size2D(2) = nlevgrnd

       isLayoutLeft = UrbanKokkosIsLayoutLeft()

       if (isLayoutLeft) then
         ! LayoutLeft: landunit index varies fastest
         idx = 0
         do j = 1, nlevgrnd
           do fc = 1, num_urbanc
             c = filter_urbanc(fc)
             if (col_pp%itype(c) == icol_road_imperv) then
               idx = idx + 1
               h2oLiq(idx) = real(h2osoi_liq(c, j), c_double)
               h2oIce(idx) = real(h2osoi_ice(c, j), c_double)
             end if
           end do
         end do
       else
         ! LayoutRight: layer index varies fastest
         idx = 0
         do fc = 1, num_urbanc
           c = filter_urbanc(fc)
           if (col_pp%itype(c) == icol_road_imperv) then
             do j = 1, nlevgrnd
               idx = idx + 1
               h2oLiq(idx) = real(h2osoi_liq(c, j), c_double)
               h2oIce(idx) = real(h2osoi_ice(c, j), c_double)
             end do
           end if
         end do
       end if

       call UrbanSetSoilLiquidWaterForImperviousRoad(urban, c_loc(h2oLiq), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetSoilIceContentForImperviousRoad(urban, c_loc(h2oIce), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       deallocate(h2oLiq)
       deallocate(h2oIce)

     end associate

   end subroutine SetImperviousRoadSoilWater

   !-----------------------------------------------------------------------
   subroutine SetRoofAndImperviousRoadTopWater(urban, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)
     !
     ! !DESCRIPTION:
     ! Seed URBANxx's TopH2OSoiLiq and TopH2OSoiIce (top layer only) for roof
     ! and impervious road from ELM's h2osoi_liq(c,1)/h2osoi_ice(c,1).
     ! Called once at startup/restart so that fwet in ComputeSurfaceFluxes
     ! starts from the correct water state rather than zero.
     !
     use ColumnType     , only : col_pp
     use ColumnDataType , only : col_ws
     use column_varcon  , only : icol_roof, icol_road_imperv
     !
     implicit none
     !
     type(UrbanType)  , intent(in) :: urban
     integer(c_int)   , intent(in) :: num_urbanl
     integer          , intent(in) :: filter_urbanl(:)
     integer(c_int)   , intent(in) :: num_urbanc
     integer          , intent(in) :: filter_urbanc(:)
     !
     integer(c_int)                       :: status
     integer                              :: fc, c
     real(c_double), allocatable, target  :: liqRoof(:), iceRoof(:)
     real(c_double), allocatable, target  :: liqImperv(:), iceImperv(:)

     associate( &
         h2osoi_liq => col_ws%h2osoi_liq, &
         h2osoi_ice => col_ws%h2osoi_ice  &
         )

       allocate(liqRoof(num_urbanl))   ; liqRoof   = 0._c_double
       allocate(iceRoof(num_urbanl))   ; iceRoof   = 0._c_double
       allocate(liqImperv(num_urbanl)) ; liqImperv = 0._c_double
       allocate(iceImperv(num_urbanl)) ; iceImperv = 0._c_double

       ! Pack top-layer values indexed by landunit (same order as filter_urbanl)
       ! Each landunit has exactly one roof column and one impervious road column.
       block
         integer, allocatable :: roofIdx(:), impIdx(:)
         integer :: fl, roofCount, impCount
         allocate(roofIdx(num_urbanl)) ; roofIdx = 0
         allocate(impIdx(num_urbanl))  ; impIdx  = 0
         roofCount = 0
         impCount  = 0
         do fc = 1, num_urbanc
           c = filter_urbanc(fc)
           select case (col_pp%itype(c))
           case (icol_roof)
             roofCount = roofCount + 1
             liqRoof(roofCount) = real(h2osoi_liq(c, 1), c_double)
             iceRoof(roofCount) = real(h2osoi_ice(c, 1), c_double)
           case (icol_road_imperv)
             impCount = impCount + 1
             liqImperv(impCount) = real(h2osoi_liq(c, 1), c_double)
             iceImperv(impCount) = real(h2osoi_ice(c, 1), c_double)
           end select
         end do
         deallocate(roofIdx, impIdx)
       end block

       call UrbanSetTopH2OSoiLiqRoof(urban, c_loc(liqRoof), num_urbanl, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetTopH2OSoiIceRoof(urban, c_loc(iceRoof), num_urbanl, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetTopH2OSoiLiqImperviousRoad(urban, c_loc(liqImperv), num_urbanl, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetTopH2OSoiIceImperviousRoad(urban, c_loc(iceImperv), num_urbanl, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       deallocate(liqRoof, iceRoof, liqImperv, iceImperv)

     end associate

     if (masterproc) then
       write(iulog,*) 'Seeded TopH2OSoiLiq/Ice for roof and impervious road from ELM data'
     end if

   end subroutine SetRoofAndImperviousRoadTopWater

   !-----------------------------------------------------------------------
   subroutine SetPerviousRoadSoilWater(urban, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)
     !
     ! !DESCRIPTION:
     ! Initialize URBANxx's persistent WaterLiquid/WaterIce for pervious road
     ! soil layers from ELM's h2osoi_liq/h2osoi_ice.  Called once at
     ! startup/restart so that heat diffusion and hydrology use correct soil
     ! moisture on the first timestep.
     !
     use elm_varpar     , only : nlevgrnd
     use ColumnType     , only : col_pp
     use ColumnDataType , only : col_ws
     use column_varcon  , only : icol_road_perv
     use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
     !
     implicit none
     !
     type(UrbanType)  , intent(in) :: urban
     integer(c_int)   , intent(in) :: num_urbanl
     integer          , intent(in) :: filter_urbanl(:)
     integer(c_int)   , intent(in) :: num_urbanc
     integer          , intent(in) :: filter_urbanc(:)
     !
     integer(c_int)                       :: status
     integer                              :: fc, c, j, idx
     integer(c_int), dimension(2)         :: size2D
     logical(c_bool)                      :: isLayoutLeft
     real(c_double), allocatable, target  :: h2oLiq(:)
     real(c_double), allocatable, target  :: h2oIce(:)

     associate( &
         h2osoi_liq => col_ws%h2osoi_liq, & ! Input: [real(r8) (:,:) ] liquid water (kg/m2)
         h2osoi_ice => col_ws%h2osoi_ice  & ! Input: [real(r8) (:,:) ] ice lens (kg/m2)
         )

       allocate(h2oLiq(num_urbanl * nlevgrnd))
       allocate(h2oIce(num_urbanl * nlevgrnd))
       size2D(1) = num_urbanl
       size2D(2) = nlevgrnd

       isLayoutLeft = UrbanKokkosIsLayoutLeft()

       if (isLayoutLeft) then
         ! LayoutLeft: landunit index varies fastest
         idx = 0
         do j = 1, nlevgrnd
           do fc = 1, num_urbanc
             c = filter_urbanc(fc)
             if (col_pp%itype(c) == icol_road_perv) then
               idx = idx + 1
               h2oLiq(idx) = real(h2osoi_liq(c, j), c_double)
               h2oIce(idx) = real(h2osoi_ice(c, j), c_double)
             end if
           end do
         end do
       else
         ! LayoutRight: layer index varies fastest
         idx = 0
         do fc = 1, num_urbanc
           c = filter_urbanc(fc)
           if (col_pp%itype(c) == icol_road_perv) then
             do j = 1, nlevgrnd
               idx = idx + 1
               h2oLiq(idx) = real(h2osoi_liq(c, j), c_double)
               h2oIce(idx) = real(h2osoi_ice(c, j), c_double)
             end do
           end if
         end do
       end if

       call UrbanSetSoilLiquidWaterForPerviousRoad(urban, c_loc(h2oLiq), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetSoilIceContentForPerviousRoad(urban, c_loc(h2oIce), size2D, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       deallocate(h2oLiq)
       deallocate(h2oIce)

     end associate

   end subroutine SetPerviousRoadSoilWater

   !-----------------------------------------------------------------------
   subroutine SetWaterTableFromState(urban, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, &
        soilhydrology_vars)
     !
     ! !DESCRIPTION:
     ! Seed URBANxx's pervious road water table depth (Zwt) and unconfined
     ! aquifer water (Wa) from ELM's soilhydrology_vars state.  Called once at
     ! startup/restart after UrbanSetup (which resets Zwt to 4.8 m and Wa to
     ! default), so that URBANxx hydrology is consistent with ELM on the first
     ! post-restart timestep.
     !
     use ColumnType        , only : col_pp
     use column_varcon     , only : icol_road_perv
     use SoilHydrologyType , only : soilhydrology_type
     !
     implicit none
     !
     type(UrbanType)         , intent(in) :: urban
     integer(c_int)          , intent(in) :: num_urbanl
     integer                 , intent(in) :: filter_urbanl(:)
     integer(c_int)          , intent(in) :: num_urbanc
     integer                 , intent(in) :: filter_urbanc(:)
     type(soilhydrology_type), intent(in) :: soilhydrology_vars
     !
     integer(c_int)                       :: status
     integer                              :: fc, c, idx
     real(c_double), allocatable, target  :: zwt_buf(:)
     real(c_double), allocatable, target  :: wa_buf(:)

     associate( &
         zwt_col => soilhydrology_vars%zwt_col , &  ! water table depth (m)
         wa_col  => soilhydrology_vars%wa_col    &  ! water in unconfined aquifer (mm)
         )

       allocate(zwt_buf(num_urbanl))
       allocate(wa_buf(num_urbanl))
       zwt_buf = 0._c_double
       wa_buf  = 0._c_double

       idx = 0
       do fc = 1, num_urbanc
         c = filter_urbanc(fc)
         if (col_pp%itype(c) == icol_road_perv) then
           idx = idx + 1
           zwt_buf(idx) = real(zwt_col(c), c_double)
           wa_buf(idx)  = real(wa_col(c),  c_double)
         end if
       end do

       call UrbanSetWaterTableDepth(urban, c_loc(zwt_buf), num_urbanl, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
       call UrbanSetAquiferWaterForPerviousRoad(urban, c_loc(wa_buf), num_urbanl, status)
       if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

       deallocate(zwt_buf)
       deallocate(wa_buf)

     end associate

   end subroutine SetWaterTableFromState

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
        urbanparams_vars, top_as, soilstate_vars, soilhydrology_vars)
     !
     implicit none
     !
     type(UrbanType)          , intent(inout) :: urban
     integer(c_int)           , intent(in)    :: num_urbanl
     integer                  , intent(in)    :: filter_urbanl(:) ! urban landunit filter
     integer(c_int)           , intent(in)    :: num_urbanc
     integer                  , intent(in)    :: filter_urbanc(:) ! urban column filter
     type(urbanparams_type)   , intent(in)    :: urbanparams_vars
     type(topounit_atmospheric_state), intent(in)    :: top_as
     type(soilstate_type)     , intent(in)    :: soilstate_vars
     type(soilhydrology_type) , intent(in)    :: soilhydrology_vars

     call SetCanyonHwr(urban, num_urbanl, filter_urbanl)
     call SetFracPervRoadOfTotalRoad(urban, num_urbanl, filter_urbanl)
     call SetWtRoof(urban, num_urbanl, filter_urbanl)
     call SetHeightParameters(urban, num_urbanl, filter_urbanl, &
          urbanparams_vars, top_as)
     call SetAlbedo(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetEmissivity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetNumberOfActiveLayersImperviousRoad(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetBuildingTemperature(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetLayerTemperatures(urban, num_urbanl, filter_urbanl)
     call SetImperviousRoadSoilWater(urban, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)
     call SetRoofAndImperviousRoadTopWater(urban, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)
     call SetPerviousRoadSoilWater(urban, num_urbanl, filter_urbanl, num_urbanc, filter_urbanc)
     call SetThermalConductivity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetHeatCapacity(urban, num_urbanl, filter_urbanl, urbanparams_vars)
     call SetSoilProperties(urban, num_urbanl, num_urbanc, filter_urbanc, soilstate_vars)
     call SetHkDepthAndTopoSlope(urban, num_urbanl, num_urbanc, filter_urbanc, soilhydrology_vars)
     call SetNMelt(urban, num_urbanl, num_urbanc, filter_urbanc)
     call SetCanyonAirStates(urban, num_urbanl, filter_urbanl)

   end subroutine SetUrbanParameters

   !-----------------------------------------------------------------------
   subroutine SetHkDepthAndTopoSlope(urban, num_urbanl, num_urbanc, filter_urbanc, soilhydrology_vars)
     !
     ! !DESCRIPTION:
     ! Set e-folding depth for saturated hydraulic conductivity and
     ! topographic slope for pervious road columns.
     ! These are time-invariant parameters set once during initialization.
     !
     use ColumnType    , only : col_pp
     use column_varcon , only : icol_road_perv
     use SoilHydrologyType, only : soilhydrology_type
     !
     implicit none
     !
     type(UrbanType)          , intent(in) :: urban
     integer(c_int)           , intent(in) :: num_urbanl
     integer(c_int)           , intent(in) :: num_urbanc
     integer                  , intent(in) :: filter_urbanc(:)
     type(soilhydrology_type) , intent(in) :: soilhydrology_vars
     !
     integer(c_int)                       :: status
     integer                              :: fc, c, idx_perv
     real(c_double), allocatable, target  :: hkdepth_buf(:)
     real(c_double), allocatable, target  :: topo_slope_buf(:)

     allocate(hkdepth_buf(num_urbanl))
     allocate(topo_slope_buf(num_urbanl))

     associate( &
          hkdepth    => soilhydrology_vars%hkdepth_col , &
          topo_slope => col_pp%topo_slope                &
          )

       idx_perv = 0
       do fc = 1, num_urbanc
         c = filter_urbanc(fc)
         if (col_pp%itype(c) == icol_road_perv) then
           idx_perv = idx_perv + 1
           hkdepth_buf(idx_perv)    = hkdepth(c)
           topo_slope_buf(idx_perv) = topo_slope(c)
         end if
       end do

     end associate

     call UrbanSetHkDepthForPerviousRoad(urban, c_loc(hkdepth_buf), num_urbanl, status)
     if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

     call UrbanSetTopoSlopeForPerviousRoad(urban, c_loc(topo_slope_buf), num_urbanl, status)
     if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

     deallocate(hkdepth_buf)
     deallocate(topo_slope_buf)

   end subroutine SetHkDepthAndTopoSlope

   !-----------------------------------------------------------------------
   subroutine SetNMelt(urban, num_urbanl, num_urbanc, filter_urbanc)
     !
     ! !DESCRIPTION:
     ! Send col_pp%n_melt (SCA shape parameter) from ELM to URBANxx for roof,
     ! impervious road, and pervious road surfaces.  Each landunit has exactly
     ! one column of each type; n_melt is a time-invariant parameter set once
     ! during initialization in initVerticalMod.
     !
     use ColumnType    , only : col_pp
     use column_varcon , only : icol_roof, icol_road_imperv, icol_road_perv
     !
     implicit none
     !
     type(UrbanType)  , intent(in) :: urban
     integer(c_int)   , intent(in) :: num_urbanl
     integer(c_int)   , intent(in) :: num_urbanc
     integer          , intent(in) :: filter_urbanc(:)
     !
     integer(c_int)                       :: status
     integer                              :: fc, c, idx_roof, idx_imperv, idx_perv
     real(c_double), allocatable, target  :: nmelt_roof(:)
     real(c_double), allocatable, target  :: nmelt_imperv(:)
     real(c_double), allocatable, target  :: nmelt_perv(:)

     allocate(nmelt_roof(num_urbanl))
     allocate(nmelt_imperv(num_urbanl))
     allocate(nmelt_perv(num_urbanl))
     nmelt_roof   = 0._c_double
     nmelt_imperv = 0._c_double
     nmelt_perv   = 0._c_double

     idx_roof   = 0
     idx_imperv = 0
     idx_perv   = 0
     do fc = 1, num_urbanc
       c = filter_urbanc(fc)
       select case (col_pp%itype(c))
       case (icol_roof)
         idx_roof = idx_roof + 1
         nmelt_roof(idx_roof) = real(col_pp%n_melt(c), c_double)
       case (icol_road_imperv)
         idx_imperv = idx_imperv + 1
         nmelt_imperv(idx_imperv) = real(col_pp%n_melt(c), c_double)
       case (icol_road_perv)
         idx_perv = idx_perv + 1
         nmelt_perv(idx_perv) = real(col_pp%n_melt(c), c_double)
       end select
     end do

     call UrbanSetNMeltRoof(urban, c_loc(nmelt_roof), num_urbanl, status)
     if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
     call UrbanSetNMeltImperviousRoad(urban, c_loc(nmelt_imperv), num_urbanl, status)
     if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
     call UrbanSetNMeltPerviousRoad(urban, c_loc(nmelt_perv), num_urbanl, status)
     if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

     deallocate(nmelt_roof)
     deallocate(nmelt_imperv)
     deallocate(nmelt_perv)

   end subroutine SetNMelt

end module UrbanxxMod
