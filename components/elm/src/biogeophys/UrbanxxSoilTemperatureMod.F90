module UrbanxxSoilTemperatureMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Set building temperature and compute heat diffusion
  ! for urban areas in the Urban++ model.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use LandunitDataType     , only : lun_es
  use UrbanxxInstanceMod   , only : urbanxx

  implicit none

  private

  ! Persistent input buffers (allocated once in init)
  real(c_double) , allocatable, target :: buildingTemp(:)
  real(c_double) , allocatable, target :: h2oLiq(:)
  real(c_double) , allocatable, target :: h2oIce(:)
  real(c_double) , allocatable, target :: h2oVol(:)

  ! Persistent output buffers for layer temperatures (allocated once in init)
  ! Roof, sunlit wall, shaded wall: (num_urbanl * nlevurb) — 1D flat arrays
  ! Impervious road, pervious road: (num_urbanl * nlevgrnd) — 1D flat arrays
  real(c_double) , allocatable, target, public :: layertemp_roof(:)
  real(c_double) , allocatable, target, public :: layertemp_improad(:)
  real(c_double) , allocatable, target, public :: layertemp_perroad(:)
  real(c_double) , allocatable, target, public :: layertemp_sunwall(:)
  real(c_double) , allocatable, target, public :: layertemp_shadwall(:)

  public :: urbanxx_soilTemperature_init
  public :: urbanxx_soilTemperature
  public :: urbanxx_soilTemperature_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilTemperature_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for heat diffusion computation.
    ! Called once during initialization.
    !
    use elm_varpar, only : nlevgrnd, nlevurb
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int) :: totalSize

    totalSize = num_urbanl * nlevgrnd

    ! Input buffers
    allocate(buildingTemp(num_urbanl))
    allocate(h2oLiq(totalSize))
    allocate(h2oIce(totalSize))
    allocate(h2oVol(totalSize))

    ! Output buffers — 2D flattened
    allocate(layertemp_roof(num_urbanl * nlevurb))
    allocate(layertemp_sunwall(num_urbanl * nlevurb))
    allocate(layertemp_shadwall(num_urbanl * nlevurb))
    allocate(layertemp_improad(num_urbanl * nlevgrnd))
    allocate(layertemp_perroad(num_urbanl * nlevgrnd))

  end subroutine urbanxx_soilTemperature_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilTemperature(num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Set building temperature and compute heat diffusion for urban areas
    !
    use TemperatureType, only : temperature_type
    use ColumnType     , only : col_pp
    use column_varcon  , only : icol_road_perv
    use elm_varpar     , only : nlevgrnd, nlevurb
    use ColumnDataType , only : col_ws
    use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int)         , intent(in) :: num_urbanl        ! number of urban landunits
    integer(c_int)         , intent(in) :: num_urbanc        ! number of urban columns
    integer                , intent(in) :: filter_urbanl(:)  ! urban layer filter
    integer                , intent(in) :: filter_urbanc(:)  ! urban column filter
    type(temperature_type) , intent(in) :: temperature_vars
    !
    ! !LOCAL VARIABLES:
    integer(c_int)                       :: status
    integer                              :: fc, c, j, fl, l
    integer                              :: idx, idx_perv, nlevbed
    integer(c_int), dimension(2)         :: size2D_urban, size2D_soil
    logical(c_bool)                      :: isLayoutLeft

    ! Set building temperature before heat diffusion
    associate(                          &
        t_building => lun_es%t_building , & ! Input: [real(r8) (:)   ]  internal building temperature (K)
        nlev2bed   => col_pp%nlevbed    , & ! Input: [integer  (:)   ]  number of layers to bedrock
        h2osoi_liq => col_ws%h2osoi_liq , & ! Input: [real(r8) (:,:)]  liquid water (kg/m2)
        h2osoi_ice => col_ws%h2osoi_ice , & ! Input: [real(r8) (:,:)]  ice lens (kg/m2)
        h2osoi_vol => col_ws%h2osoi_vol   & ! Input: [real(r8) (:,:)]  volumetric soil water [m3/m3]
        )

      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         buildingTemp(fl) = t_building(l)
      end do

      call UrbanSetBuildingTemperature(urbanxx, c_loc(buildingTemp), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! Set soil water content for pervious road so heat diffusion uses
      ! up-to-date water state (thermal conductivity depends on water content)
      size2D_soil(1) = num_urbanl
      size2D_soil(2) = nlevgrnd

      isLayoutLeft = UrbanKokkosIsLayoutLeft()

      if (isLayoutLeft) then
        ! LayoutLeft: First dimension (landunits) varies fastest
        idx = 0
        do j = 1, nlevgrnd
          idx_perv = 0
          do fc = 1, num_urbanc
            c = filter_urbanc(fc)
            if (col_pp%itype(c) == icol_road_perv) then
              idx_perv = idx_perv + 1
              idx = idx + 1
              nlevbed = nlev2bed(c)
              if (j <= nlevbed) then
                h2oLiq(idx) = h2osoi_liq(c, j)
                h2oIce(idx) = h2osoi_ice(c, j)
                h2oVol(idx) = h2osoi_vol(c, j)
              else
                h2oLiq(idx) = 0.0_r8
                h2oIce(idx) = 0.0_r8
                h2oVol(idx) = 0.0_r8
              end if
            end if
          end do
        end do
      else
        ! LayoutRight: Last dimension (layers) varies fastest
        idx = 0
        do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col_pp%itype(c) == icol_road_perv) then
            nlevbed = nlev2bed(c)
            do j = 1, nlevgrnd
              idx = idx + 1
              if (j <= nlevbed) then
                h2oLiq(idx) = h2osoi_liq(c, j)
                h2oIce(idx) = h2osoi_ice(c, j)
                h2oVol(idx) = h2osoi_vol(c, j)
              else
                h2oLiq(idx) = 0.0_r8
                h2oIce(idx) = 0.0_r8
                h2oVol(idx) = 0.0_r8
              end if
            end do
          end if
        end do
      end if

      call UrbanSetSoilLiquidWaterForPerviousRoad(urbanxx, c_loc(h2oLiq), size2D_soil, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetSoilIceContentForPerviousRoad(urbanxx, c_loc(h2oIce), size2D_soil, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetSoilVolumetricWaterForPerviousRoad(urbanxx, c_loc(h2oVol), size2D_soil, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanComputeHeatDiffusion(urbanxx, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! Extract layer temperatures from UrbanXX
      ! Roof, sunlit wall, shaded wall: (num_urbanl, nlevurb)
      size2D_urban(1) = num_urbanl
      size2D_urban(2) = nlevurb

      call UrbanGetLayerTempRoof(urbanxx, c_loc(layertemp_roof), size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetLayerTempSunlitWall(urbanxx, c_loc(layertemp_sunwall), size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetLayerTempShadedWall(urbanxx, c_loc(layertemp_shadwall), size2D_urban, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! Impervious road, pervious road: (num_urbanl, nlevgrnd)
      size2D_soil(1) = num_urbanl
      size2D_soil(2) = nlevgrnd

      call UrbanGetLayerTempImperviousRoad(urbanxx, c_loc(layertemp_improad), size2D_soil, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetLayerTempPerviousRoad(urbanxx, c_loc(layertemp_perroad), size2D_soil, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    end associate

  end subroutine urbanxx_soilTemperature

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilTemperature_check(num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Set building temperature and compute heat diffusion for urban areas
    !
    use TemperatureType, only : temperature_type
    use ColumnType     , only : col_pp
    use column_varcon  , only : icol_road_perv
    use elm_varpar     , only : nlevgrnd, nlevurb
    use column_varcon  , only : icol_roof, icol_road_imperv, icol_road_perv, icol_sunwall, icol_shadewall
    use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
    use ColumnDataType    , only : col_es
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int)         , intent(in) :: num_urbanl        ! number of urban landunits
    integer(c_int)         , intent(in) :: num_urbanc        ! number of urban columns
    integer                , intent(in) :: filter_urbanl(:)  ! urban layer filter
    integer                , intent(in) :: filter_urbanc(:)  ! urban column filter
    type(temperature_type) , intent(in) :: temperature_vars
    real(c_double) , allocatable, target :: layertemp_roof_2d(:,:)
    real(c_double) , allocatable, target :: layertemp_sunwall_2d(:,:)
    real(c_double) , allocatable, target :: layertemp_shadwall_2d(:,:)
    real(c_double) , allocatable, target :: layertemp_improad_2d(:,:)
    real(c_double) , allocatable, target :: layertemp_perroad_2d(:,:)

    !
    ! !LOCAL VARIABLES:
    integer(c_int)                       :: status
    integer                              :: fc, c, j, fl, l
    integer(c_int), dimension(2)         :: size2D_urban, size2D_soil
    real(r8)                             :: max_error, max_rel_error
    integer :: idx_roof, idx_road_imperv, idx_road_perv, idx_sunwall, idx_shadwall, idx_landunit
    integer :: count

    associate(                          &
         t_soisno                => col_es%t_soisno  & ! Output: [real(r8) (:,:) ]  soil temperature (Kelvin)
        )
    max_error     = 0._r8
    max_rel_error = 0._r8
    idx_roof = 0
    idx_road_imperv = 0
    idx_road_perv = 0
    idx_sunwall = 0
    idx_shadwall = 0

      allocate(layertemp_roof_2d(num_urbanl, nlevurb))
      allocate(layertemp_sunwall_2d(num_urbanl, nlevurb))
      allocate(layertemp_shadwall_2d(num_urbanl, nlevurb))
      allocate(layertemp_improad_2d(num_urbanl, nlevgrnd))
      allocate(layertemp_perroad_2d(num_urbanl, nlevgrnd))

      if (UrbanKokkosIsLayoutLeft()) then
        count = 0
        do j = 1, nlevurb
          do l = 1, num_urbanl
            count = count + 1
            layertemp_roof_2d(l,j) = layertemp_roof(count)
            layertemp_sunwall_2d(l,j) = layertemp_sunwall(count)
            layertemp_shadwall_2d(l,j) = layertemp_shadwall(count)
          end do
        end do
        count = 0
        do j = 1, nlevgrnd
          do l = 1, num_urbanl
            count = count + 1
            layertemp_improad_2d(l,j) = layertemp_improad(count)
            layertemp_perroad_2d(l,j) = layertemp_perroad(count)
          end do
        end do
      else
        count = 0
        do l = 1, num_urbanl
          do j = 1, nlevurb
            count = count + 1
            layertemp_roof_2d(l,j) = layertemp_roof(count)
            layertemp_sunwall_2d(l,j) = layertemp_sunwall(count)
            layertemp_shadwall_2d(l,j) = layertemp_shadwall(count)
          end do
        end do
        count = 0
        do l = 1, num_urbanl
          do j = 1, nlevgrnd
            count = count + 1
            layertemp_improad_2d(l,j) = layertemp_improad(count)
            layertemp_perroad_2d(l,j) = layertemp_perroad(count)
          end do
        end do
      end if

      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        select case (col_pp%itype(c))
        case (icol_roof)
          idx_roof = idx_roof + 1
          do j = 1, nlevurb
            max_error     = max(max_error,     abs(t_soisno(c,j) - layertemp_roof_2d(idx_roof,j)))
            max_rel_error = max(max_rel_error, abs(t_soisno(c,j) - layertemp_roof_2d(idx_roof,j))     / max(abs(t_soisno(c,j)), 1.0e-20_r8))
          end do
        case (icol_sunwall)
          idx_sunwall = idx_sunwall + 1
          do j = 1, nlevurb
            max_error     = max(max_error,     abs(t_soisno(c,j) - layertemp_sunwall_2d(idx_sunwall,j)))
            max_rel_error = max(max_rel_error, abs(t_soisno(c,j) - layertemp_sunwall_2d(idx_sunwall,j)) / max(abs(t_soisno(c,j)), 1.0e-20_r8))
          end do
        case (icol_shadewall)
          idx_shadwall = idx_shadwall + 1
          do j = 1, nlevurb
            max_error     = max(max_error,     abs(t_soisno(c,j) - layertemp_shadwall_2d(idx_shadwall,j)))
            max_rel_error = max(max_rel_error, abs(t_soisno(c,j) - layertemp_shadwall_2d(idx_shadwall,j)) / max(abs(t_soisno(c,j)), 1.0e-20_r8))
          end do
        case (icol_road_imperv)
          idx_road_imperv = idx_road_imperv + 1
          do j = 1, nlevgrnd
            max_error     = max(max_error,     abs(t_soisno(c,j) - layertemp_improad_2d(idx_road_imperv,j)))
            max_rel_error = max(max_rel_error, abs(t_soisno(c,j) - layertemp_improad_2d(idx_road_imperv,j)) / max(abs(t_soisno(c,j)), 1.0e-20_r8))
          end do
        case (icol_road_perv)
          idx_road_perv = idx_road_perv + 1
          do j = 1, nlevgrnd
            max_error     = max(max_error,     abs(t_soisno(c,j) - layertemp_perroad_2d(idx_road_perv,j)))
            max_rel_error = max(max_rel_error, abs(t_soisno(c,j) - layertemp_perroad_2d(idx_road_perv,j)) / max(abs(t_soisno(c,j)), 1.0e-20_r8))
          end do
        end select
        if (max_error > 1.0e-9) then
          write(iulog,*)'Max error in temperature: ', max_error, ' at column ', c, ' itype = ', col_pp%itype(c), &
          'exceed tolerance'
          call exit(0)
        endif
      end do
      write(iulog,*) 'Max error in soil temperatures      : ', max_error, ' (rel: ', max_rel_error, ')'

      deallocate(layertemp_roof_2d)
      deallocate(layertemp_sunwall_2d)
      deallocate(layertemp_shadwall_2d)
      deallocate(layertemp_improad_2d)
      deallocate(layertemp_perroad_2d)

    end associate
  end subroutine urbanxx_soilTemperature_check

  end module UrbanxxSoilTemperatureMod
