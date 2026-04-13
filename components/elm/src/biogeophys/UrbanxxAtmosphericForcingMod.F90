module UrbanxxAtmosphericForcingMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Set atmospheric forcing data for the Urban++ model.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use urban_kokkos_interface
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use UrbanParamsType      , only : urbanparams_type
  use SurfaceAlbedoType    , only : surfalb_type
  use LandunitType         , only : lun_pp
  use FrictionVelocityType , only : frictionvel_type
  use TopounitDataType     , only : top_as, top_af
  use UrbanxxInstanceMod   , only : urbanxx, numBands, numTypes
  use UrbanxxMod           , only : SetHeightParameters
  use TopounitDataType     , only : topounit_atmospheric_state
  use GridcellType         , only : grc_pp

  implicit none

  private

  ! Persistent input buffers (allocated once in init)
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

  public :: urbanxx_SetAtmosphericForcing_init
  public :: urbanxx_SetAtmosphericForcing

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_SetAtmosphericForcing_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for atmospheric forcing.
    ! Called once during initialization.
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int) :: totalSize3D

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

    totalSize3D = num_urbanl * numBands * numTypes
    allocate(atmShortwave(totalSize3D))

  end subroutine urbanxx_SetAtmosphericForcing_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_SetAtmosphericForcing(num_urbanl, filter_urbanl, nextsw_cday  , declinp1, surfalb_vars, &
         urbanparams_vars, top_as)
    !
    use shr_orb_mod, only : shr_orb_cosz
    !
    implicit none
    !
    integer(c_int)     , intent(in) :: num_urbanl
    integer            , intent(in) :: filter_urbanl(:)         ! urban landunit filter
    type(surfalb_type) , intent(in) :: surfalb_vars
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(topounit_atmospheric_state) , intent(in)    :: top_as

    !
    integer(c_int)                       :: status
    integer                              :: fl, l, t, g, iband, itype, idx
    integer(c_int)                       :: totalSize3D
    integer(c_int), dimension(3)         :: size3D
    logical(c_bool)                      :: isLayoutLeft
    real(r8), intent(in) :: nextsw_cday        ! calendar day at Greenwich (1.00, ..., days/year)
    real(r8), intent(in) :: declinp1           ! declination angle (radians) for next time step

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

      size3D = [num_urbanl, numBands, numTypes]
      totalSize3D = num_urbanl * numBands * numTypes

      ! Fill arrays with values from ELM data structures
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         t = lun_pp%topounit(l)
         g = lun_pp%gridcell(l)

         atmTemp(fl)     = forc_t(t)
         atmPotTemp(fl)  = forc_th(t)
         atmRho(fl)      = forc_rho(t)
         atmSpcHumd(fl)  = forc_q(t)
         atmPress(fl)    = forc_pbot(t)
         atmWindU(fl)    = forc_u(t)
         atmWindV(fl)    = forc_v(t)
         atmCoszen(fl)   = shr_orb_cosz (nextsw_cday, grc_pp%lat(g), grc_pp%lon(g), declinp1)
         atmFracSnow(fl) = forc_snow(t)
         atmLongwave(fl) = forc_lwrad(t)
      end do

      ! Fill shortwave arrays with direct and diffuse for VIS and NIR bands
      ! Kokkos View dimensions: (numLandunits, numRadBands, numRadTypes)
      ! itype = 0: diffuse, itype = 1: direct
      isLayoutLeft = UrbanKokkosIsLayoutLeft()

      if (isLayoutLeft) then
        ! LayoutLeft: First dimension (landunits) varies fastest
        ! Memory order: landunit, band, type
        idx = 0
        do itype = 0, numTypes - 1
          do iband = 0, numBands - 1
            do fl = 1, num_urbanl
              l = filter_urbanl(fl)
              t = lun_pp%topounit(l)
              idx = idx + 1
              if (itype == 0) then
                atmShortwave(idx) = forc_solad(t, iband+1)  ! direct
              else
                atmShortwave(idx) = forc_solai(t, iband+1)  ! diffuse
              end if
            end do
          end do
        end do
      else
        ! LayoutRight: Last dimension (types) varies fastest
        ! Memory order: type, band, landunit
        idx = 0
        do fl = 1, num_urbanl
          l = filter_urbanl(fl)
          t = lun_pp%topounit(l)
          do iband = 0, numBands - 1
            do itype = 0, numTypes - 1
              idx = idx + 1
              if (itype == 0) then
                atmShortwave(idx) = forc_solad(t, iband+1)  ! direct
              else
                atmShortwave(idx) = forc_solai(t, iband+1)  ! diffuse
              end if
            end do
          end do
        end do
      end if

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

      call UrbanComputeNetShortwaveRadiation(urbanxx, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      if (masterproc) then
         write(iulog,*) 'Set atmospheric forcing from ELM data structures'
      end if

      call SetHeightParameters(urbanxx, num_urbanl, filter_urbanl, &
         urbanparams_vars, top_as)

    end associate

  end subroutine urbanxx_SetAtmosphericForcing

end module UrbanxxAtmosphericForcingMod
