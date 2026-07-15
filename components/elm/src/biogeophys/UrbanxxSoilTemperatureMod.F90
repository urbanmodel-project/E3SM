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
  use UrbanxxInstanceMod   , only : urbanxx

  implicit none

  private

  ! Persistent input buffers (allocated once in init)
  ! Persistent output buffers for layer temperatures (allocated once in init)
  ! Roof, sunlit wall, shaded wall: (num_urbanl * nlevurb) — 1D flat arrays
  ! Impervious road, pervious road: (num_urbanl * nlevgrnd) — 1D flat arrays
  real(c_double) , allocatable, target, public :: layertemp_roof(:)
  real(c_double) , allocatable, target, public :: layertemp_improad(:)
  real(c_double) , allocatable, target, public :: layertemp_perroad(:)
  real(c_double) , allocatable, target, public :: layertemp_sunwall(:)
  real(c_double) , allocatable, target, public :: layertemp_shadwall(:)

  ! Native URBANxx snow-melt flux, for validation against ELM's own
  ! col_wf%qflx_snomelt in urbanxx_soilTemperature_check — (num_urbanl) 1D
  real(c_double) , allocatable, target, public :: qflxsnowmelt_roof(:)
  real(c_double) , allocatable, target, public :: qflxsnowmelt_improad(:)
  real(c_double) , allocatable, target, public :: qflxsnowmelt_perroad(:)

  public :: urbanxx_soilTemperature_init
  public :: urbanxx_soilTemperature
  public :: urbanxx_soilTemperature_check

  logical, save :: mapping_reported = .false.

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

    ! Output buffers — 2D flattened
    allocate(layertemp_roof(num_urbanl * nlevurb))
    allocate(layertemp_sunwall(num_urbanl * nlevurb))
    allocate(layertemp_shadwall(num_urbanl * nlevurb))
    allocate(layertemp_improad(num_urbanl * nlevgrnd))
    allocate(layertemp_perroad(num_urbanl * nlevgrnd))

    allocate(qflxsnowmelt_roof(num_urbanl))
    allocate(qflxsnowmelt_improad(num_urbanl))
    allocate(qflxsnowmelt_perroad(num_urbanl))

  end subroutine urbanxx_soilTemperature_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilTemperature(num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Compute heat diffusion for urban areas.
    ! Impervious road soil moisture (WaterLiquid/WaterIce) is maintained
    ! as persistent URBANxx state, initialized once from ELM at startup
    ! and updated each timestep via an internal phase change step inside
    ! UrbanComputeHeatDiffusion.  No per-timestep data transfer from ELM.
    !
    use TemperatureType, only : temperature_type
    use ColumnType     , only : col_pp
    use ColumnDataType , only : col_es, col_ws
    use column_varcon  , only : icol_road_perv
    use elm_varpar     , only : nlevgrnd, nlevurb
    use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
    use elm_time_manager  , only : get_nstep
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
    integer                              :: fc, c, j, fl, l, count
    integer                              :: idx_road_perv
    integer(c_int), dimension(2)         :: size2D_urban, size2D_soil
    real(c_double), allocatable, target  :: layertemp_perroad_2d(:,:)
    real(c_double), allocatable, target  :: waterliq_perroad_2d(:,:)
    real(c_double), allocatable, target  :: waterice_perroad_2d(:,:)
    associate(                        &
        nlev2bed => col_pp%nlevbed    , & ! Input: [integer  (:)   ]  number of layers to bedrock
        t_soisno => col_es%t_soisno   , & ! Input: [real(r8) (:,:) ]  soil temperature (Kelvin)
        h2osoi_liq => col_ws%h2osoi_liq, & ! Input: [real(r8) (:,:) ] liquid water (kg/m2)
        h2osoi_ice => col_ws%h2osoi_ice  & ! Input: [real(r8) (:,:) ] ice water (kg/m2)
        )

      if (get_nstep() == 1) then
        size2D_soil(1) = num_urbanl
        size2D_soil(2) = nlevgrnd

        call UrbanGetLayerTempPerviousRoad(urbanxx, c_loc(layertemp_perroad), size2D_soil, status)
        if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

        allocate(layertemp_perroad_2d(num_urbanl, nlevgrnd))
        allocate(waterliq_perroad_2d(num_urbanl, nlevgrnd))
        allocate(waterice_perroad_2d(num_urbanl, nlevgrnd))

        call UrbanGetSoilLiquidWaterPerviousRoad(urbanxx, c_loc(waterliq_perroad_2d), size2D_soil, status)
        if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
        call UrbanGetSoilIceContentPerviousRoad(urbanxx, c_loc(waterice_perroad_2d), size2D_soil, status)
        if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

        if (UrbanKokkosIsLayoutLeft()) then
          count = 0
          do j = 1, nlevgrnd
            do l = 1, num_urbanl
              count = count + 1
              layertemp_perroad_2d(l,j) = layertemp_perroad(count)
            end do
          end do
        else
          count = 0
          do l = 1, num_urbanl
            do j = 1, nlevgrnd
              count = count + 1
              layertemp_perroad_2d(l,j) = layertemp_perroad(count)
            end do
          end do
        end if

        idx_road_perv = 0
        do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col_pp%itype(c) == icol_road_perv) then
            idx_road_perv = idx_road_perv + 1
            if (idx_road_perv == 1) then
              write(iulog,'(A,I0,A,I0)') 'URBANxx pre-solve pervious road state: column=', c, ', packed_index=', idx_road_perv
              do j = 1, nlevgrnd
                write(iulog,'(A,I0,6(A,ES24.16))') '  layer=', j, &
                     ', elm_t=', t_soisno(c,j), ', uxx_t=', layertemp_perroad_2d(idx_road_perv,j), &
                     ', elm_liq=', h2osoi_liq(c,j), ', uxx_liq=', waterliq_perroad_2d(idx_road_perv,j), &
                     ', elm_ice=', h2osoi_ice(c,j), ', uxx_ice=', waterice_perroad_2d(idx_road_perv,j)
              end do
              exit
            end if
          end if
        end do

        deallocate(layertemp_perroad_2d)
        deallocate(waterliq_perroad_2d)
        deallocate(waterice_perroad_2d)
      end if

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

      ! Native URBANxx snow-melt flux (already computed above, inside
      ! UrbanComputeHeatDiffusion's phase-change step) — retrieved here for
      ! validation against ELM's own qflx_snomelt in
      ! urbanxx_soilTemperature_check.
      call UrbanGetQflxSnowMeltRoof(urbanxx, c_loc(qflxsnowmelt_roof), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetQflxSnowMeltImperviousRoad(urbanxx, c_loc(qflxsnowmelt_improad), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetQflxSnowMeltPerviousRoad(urbanxx, c_loc(qflxsnowmelt_perroad), num_urbanl, status)
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
    use ColumnDataType    , only : col_es, col_ws, col_ef, col_wf
    use elm_time_manager  , only : get_nstep, get_curr_date
    use UrbanxxDiagnosticControlMod, only : urbanxx_softfail_permits
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
    real(r8)                             :: column_max_error, column_max_rel_error
    real(r8)                             :: elm_value, urbanxx_value, signed_diff, abs_diff, rel_diff
    real(r8)                             :: flux_max_error, flux_max_rel_error
    real(r8)                             :: column_flux_error, column_flux_rel_error
    real(r8)                             :: flux_elm_value, flux_urbanxx_value
    real(r8)                             :: flux_signed_diff, flux_abs_diff, flux_rel_diff
    integer :: idx_roof, idx_road_imperv, idx_road_perv, idx_sunwall, idx_shadwall, idx_landunit
    integer :: count, column_max_layer, column_packed_idx
    integer :: year, month, day, tod

    associate(                          &
         t_soisno                => col_es%t_soisno, & ! Output: [real(r8) (:,:) ]  soil temperature (Kelvin)
         h2osno                  => col_ws%h2osno, & ! Input: [real(r8) (:)] snow water (kg/m2)
         snow_depth              => col_ws%snow_depth, & ! Input: [real(r8) (:)] snow depth (m)
         frac_sno                => col_ws%frac_sno, & ! Input: [real(r8) (:)] snow fraction (-)
         h2osoi_liq              => col_ws%h2osoi_liq, & ! Input: [real(r8) (:,:)] liquid water (kg/m2)
         h2osoi_ice              => col_ws%h2osoi_ice, & ! Input: [real(r8) (:,:)] ice water (kg/m2)
         eflx_dhsdT              => col_ef%eflx_dhsdT, & ! Input: [real(r8) (:)] derivative of hs wrt temperature
         eflx_hs_top_snow        => col_ef%eflx_hs_top_snow, & ! Input: [real(r8) (:)] heat flux on top snow layer
         eflx_sabg_lyr           => col_ef%eflx_sabg_lyr, & ! Input: [real(r8) (:,:)] absorbed solar by layer
         htvp                    => col_ef%htvp, & ! Input: [real(r8) (:)] latent heat
         qflx_snomelt            => col_wf%qflx_snomelt  & ! Input: [real(r8) (:)] snow melt (mm H2O/s)
        )
    max_error     = 0._r8
    max_rel_error = 0._r8
    flux_max_error     = 0._r8
    flux_max_rel_error = 0._r8
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
        idx_landunit = col_pp%landunit(c)
        column_max_error = 0._r8
        column_max_rel_error = 0._r8
        column_max_layer = -1
        column_packed_idx = -1
        column_flux_error = 0._r8
        column_flux_rel_error = 0._r8
        select case (col_pp%itype(c))
        case (icol_roof)
          idx_roof = idx_roof + 1
          column_packed_idx = idx_roof
          do j = 1, nlevurb
            call record_error(t_soisno(c,j), layertemp_roof_2d(idx_roof,j), j)
          end do
          call record_flux_error(qflx_snomelt(c), qflxsnowmelt_roof(idx_roof))
        case (icol_sunwall)
          idx_sunwall = idx_sunwall + 1
          column_packed_idx = idx_sunwall
          do j = 1, nlevurb
            call record_error(t_soisno(c,j), layertemp_sunwall_2d(idx_sunwall,j), j)
          end do
        case (icol_shadewall)
          idx_shadwall = idx_shadwall + 1
          column_packed_idx = idx_shadwall
          do j = 1, nlevurb
            call record_error(t_soisno(c,j), layertemp_shadwall_2d(idx_shadwall,j), j)
          end do
        case (icol_road_imperv)
          idx_road_imperv = idx_road_imperv + 1
          column_packed_idx = idx_road_imperv
          do j = 1, nlevgrnd
            call record_error(t_soisno(c,j), layertemp_improad_2d(idx_road_imperv,j), j)
          end do
          call record_flux_error(qflx_snomelt(c), qflxsnowmelt_improad(idx_road_imperv))
        case (icol_road_perv)
          idx_road_perv = idx_road_perv + 1
          column_packed_idx = idx_road_perv
          if (.not. mapping_reported) then
            write(iulog,'(A,4(I0,A))') 'URBANxx pervious-road mapping: column=', c, ', landunit=', idx_landunit, &
                 ', packed_index=', idx_road_perv, ', column_type=', col_pp%itype(c), ''
          end if
          do j = 1, nlevgrnd
            call record_error(t_soisno(c,j), layertemp_perroad_2d(idx_road_perv,j), j)
          end do
          call record_flux_error(qflx_snomelt(c), qflxsnowmelt_perroad(idx_road_perv))
        end select

        if (column_max_error > 1.0e-9_r8) then
          call get_curr_date(year, month, day, tod)
          write(iulog,'(A)') 'ERROR: URBANxx soil-temperature comparison exceeds tolerance'
          write(iulog,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '  step=', get_nstep(), ', date=', year, '-', month, '-', day, ', tod=', tod
          write(iulog,'(A,I0,A,I0,A,I0,A,I0)') '  column=', c, ', landunit=', idx_landunit, &
               ', packed_index=', column_packed_idx, ', column_type=', col_pp%itype(c)
          write(iulog,'(A,I0)') '  layer=', column_max_layer
          write(iulog,'(A,ES24.16)') '  ELM temperature       =', elm_value
          write(iulog,'(A,ES24.16)') '  URBANxx temperature   =', urbanxx_value
          write(iulog,'(A,ES24.16)') '  signed difference     =', signed_diff
          write(iulog,'(A,ES24.16)') '  absolute difference   =', column_max_error
          write(iulog,'(A,ES24.16)') '  relative difference   =', column_max_rel_error
          write(iulog,'(A,3ES24.16)') '  ELM snow state        =', h2osno(c), snow_depth(c), frac_sno(c)
          write(iulog,'(A,2ES24.16)') '  ELM top liq/ice       =', h2osoi_liq(c,1), h2osoi_ice(c,1)
          write(iulog,'(A,3ES24.16)') '  ELM flux deriv/state  =', eflx_dhsdT(c), eflx_hs_top_snow(c), htvp(c)
          write(iulog,'(A,ES24.16)') '  ELM sabg_lyr top      =', eflx_sabg_lyr(c,1)
          if (urbanxx_softfail_permits(column_max_error, 1.0e-9_r8)) then
            write(iulog,'(A)') '  URBANXX_SOFT_FAIL: exceedance within ceiling; logging and continuing'
          else
            call exit(0)
          endif
        endif

        if (column_flux_error > 1.0e-9_r8) then
          call get_curr_date(year, month, day, tod)
          write(iulog,'(A)') 'ERROR: URBANxx snow-melt flux comparison exceeds tolerance'
          write(iulog,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '  step=', get_nstep(), ', date=', year, '-', month, '-', day, ', tod=', tod
          write(iulog,'(A,I0,A,I0,A,I0,A,I0)') '  column=', c, ', landunit=', idx_landunit, &
               ', packed_index=', column_packed_idx, ', column_type=', col_pp%itype(c)
          write(iulog,'(A,I0)') '  ELM snl (snow layers) =', col_pp%snl(c)
          write(iulog,'(A,ES24.16)') '  ELM qflx_snomelt      =', flux_elm_value
          write(iulog,'(A,ES24.16)') '  URBANxx QflxSnowMelt  =', flux_urbanxx_value
          write(iulog,'(A,ES24.16)') '  signed difference     =', flux_signed_diff
          write(iulog,'(A,ES24.16)') '  absolute difference   =', column_flux_error
          write(iulog,'(A,ES24.16)') '  relative difference   =', column_flux_rel_error
          write(iulog,'(A,3ES24.16)') '  ELM snow state        =', h2osno(c), snow_depth(c), frac_sno(c)
          write(iulog,'(A,2ES24.16)') '  ELM top liq/ice       =', h2osoi_liq(c,1), h2osoi_ice(c,1)
          if (urbanxx_softfail_permits(column_flux_error, 1.0e-9_r8)) then
            write(iulog,'(A)') '  URBANXX_SOFT_FAIL: exceedance within ceiling; logging and continuing'
          else
            call exit(0)
          endif
        endif
      end do
      mapping_reported = .true.
      write(iulog,*) 'Max error in soil temperatures      : ', max_error, ' (rel: ', max_rel_error, ')'
      write(iulog,*) 'Max error in snow-melt flux         : ', flux_max_error, ' (rel: ', flux_max_rel_error, ')'

      deallocate(layertemp_roof_2d)
      deallocate(layertemp_sunwall_2d)
      deallocate(layertemp_shadwall_2d)
      deallocate(layertemp_improad_2d)
      deallocate(layertemp_perroad_2d)

    end associate

  contains

    subroutine record_error(elm_temperature, urbanxx_temperature, layer)
      real(r8), intent(in) :: elm_temperature, urbanxx_temperature
      integer, intent(in) :: layer

      abs_diff = abs(elm_temperature - urbanxx_temperature)
      rel_diff = abs_diff / max(abs(elm_temperature), 1.0e-20_r8)
      max_rel_error = max(max_rel_error, rel_diff)
      if (abs_diff > column_max_error) then
        column_max_error = abs_diff
        column_max_rel_error = rel_diff
        column_max_layer = layer
        elm_value = elm_temperature
        urbanxx_value = urbanxx_temperature
        signed_diff = elm_temperature - urbanxx_temperature
      end if
      if (abs_diff > max_error) then
        max_error = abs_diff
      end if
    end subroutine record_error

    subroutine record_flux_error(elm_flux, urbanxx_flux)
      real(r8), intent(in) :: elm_flux, urbanxx_flux

      flux_abs_diff = abs(elm_flux - urbanxx_flux)
      flux_rel_diff = flux_abs_diff / max(abs(elm_flux), 1.0e-20_r8)
      flux_max_rel_error = max(flux_max_rel_error, flux_rel_diff)
      if (flux_abs_diff > column_flux_error) then
        column_flux_error = flux_abs_diff
        column_flux_rel_error = flux_rel_diff
        flux_elm_value = elm_flux
        flux_urbanxx_value = urbanxx_flux
        flux_signed_diff = elm_flux - urbanxx_flux
      end if
      if (flux_abs_diff > flux_max_error) then
        flux_max_error = flux_abs_diff
      end if
    end subroutine record_flux_error
  end subroutine urbanxx_soilTemperature_check

  end module UrbanxxSoilTemperatureMod
