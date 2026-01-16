module UrbanxxMod

  use iso_c_binding
  use urban_mod
  use urban_kokkos_interface
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use spmdMod           , only : masterproc, iam
  use shr_sys_mod       , only : shr_sys_flush
  use elm_varctl        , only : iulog
  use decompMod         , only : bounds_type
  use UrbanParamsType   , only : urbanparams_type
  use SolarAbsorbedType , only : solarabs_type
  use SurfaceAlbedoType , only : surfalb_type
  use LandunitType      , only : lun_pp

  implicit none

  private

  type(UrbanType) :: urbanxx

  public :: urbanxx_initialize

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_initialize(bounds, num_urbanl, filter_urbanl, &
       num_urbanc, filter_urbanc, num_urbanp, filter_urbanp, &
       urbanparams_vars, solarabs_vars, surfalb_vars)
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

    call SetUrbanParameters(urbanxx, num_urbanl, filter_urbanl)

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
  subroutine SetUrbanParameters(urban, num_urbanl, filter_urbanl)
    !
    implicit none
    !
    type(UrbanType), intent(inout) :: urban
    integer(c_int) , intent(in)    :: num_urbanl
    integer        , intent(in)    :: filter_urbanl(:) ! urban landunit filter

    call SetCanyonHwr(urban, num_urbanl, filter_urbanl)
    call SetFracPervRoadOfTotalRoad(urban, num_urbanl, filter_urbanl)
    call SetWtRoof(urban, num_urbanl, filter_urbanl)

  end subroutine SetUrbanParameters

end module UrbanxxMod
