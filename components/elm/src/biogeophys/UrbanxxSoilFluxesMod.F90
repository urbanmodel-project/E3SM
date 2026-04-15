module UrbanxxSoilFluxesMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute soil fluxes (ground heat flux, evaporation/dew partition) for
  ! urban areas using the Urban++ model.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use UrbanxxInstanceMod   , only : urbanxx
  use VegetationDataType   , only : veg_ef, veg_wf
  use VegetationType       , only : veg_pp

  implicit none

  private

  ! Persistent output buffers (allocated once in init)
  real(c_double) , allocatable, target :: eflx_soil_grnd_roof(:)
  real(c_double) , allocatable, target :: eflx_soil_grnd_improad(:)
  real(c_double) , allocatable, target :: eflx_soil_grnd_perroad(:)
  real(c_double) , allocatable, target :: eflx_soil_grnd_sunwall(:)
  real(c_double) , allocatable, target :: eflx_soil_grnd_shadwall(:)
  real(c_double) , allocatable, target :: qflx_evap_grnd_roof(:)
  real(c_double) , allocatable, target :: qflx_evap_grnd_improad(:)
  real(c_double) , allocatable, target :: qflx_evap_grnd_perroad(:)
  real(c_double) , allocatable, target :: qflx_sub_snow_roof(:)
  real(c_double) , allocatable, target :: qflx_sub_snow_improad(:)
  real(c_double) , allocatable, target :: qflx_sub_snow_perroad(:)
  real(c_double) , allocatable, target :: qflx_dew_snow_roof(:)
  real(c_double) , allocatable, target :: qflx_dew_snow_improad(:)
  real(c_double) , allocatable, target :: qflx_dew_snow_perroad(:)
  real(c_double) , allocatable, target :: qflx_dew_grnd_roof(:)
  real(c_double) , allocatable, target :: qflx_dew_grnd_improad(:)
  real(c_double) , allocatable, target :: qflx_dew_grnd_perroad(:)

  ! Input buffers: top-layer soil water for roof and impervious road
  real(c_double) , allocatable, target :: top_h2osoi_liq_roof(:)
  real(c_double) , allocatable, target :: top_h2osoi_ice_roof(:)
  real(c_double) , allocatable, target :: top_h2osoi_liq_improad(:)
  real(c_double) , allocatable, target :: top_h2osoi_ice_improad(:)

  public :: urbanxx_soilFluxes_init
  public :: urbanxx_soilFluxes
  public :: urbanxx_soilFluxes_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilFluxes_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for soil fluxes computation.
    ! Called once during initialization.
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl

    ! Ground heat flux — all five surfaces
    allocate(eflx_soil_grnd_roof(num_urbanl))
    allocate(eflx_soil_grnd_improad(num_urbanl))
    allocate(eflx_soil_grnd_perroad(num_urbanl))
    allocate(eflx_soil_grnd_sunwall(num_urbanl))
    allocate(eflx_soil_grnd_shadwall(num_urbanl))

    ! Evaporation/dew partition — roof, impervious road, pervious road
    allocate(qflx_evap_grnd_roof(num_urbanl))
    allocate(qflx_evap_grnd_improad(num_urbanl))
    allocate(qflx_evap_grnd_perroad(num_urbanl))
    allocate(qflx_sub_snow_roof(num_urbanl))
    allocate(qflx_sub_snow_improad(num_urbanl))
    allocate(qflx_sub_snow_perroad(num_urbanl))
    allocate(qflx_dew_snow_roof(num_urbanl))
    allocate(qflx_dew_snow_improad(num_urbanl))
    allocate(qflx_dew_snow_perroad(num_urbanl))
    allocate(qflx_dew_grnd_roof(num_urbanl))
    allocate(qflx_dew_grnd_improad(num_urbanl))
    allocate(qflx_dew_grnd_perroad(num_urbanl))

    ! Top-layer soil water inputs for roof and impervious road
    allocate(top_h2osoi_liq_roof(num_urbanl))
    allocate(top_h2osoi_ice_roof(num_urbanl))
    allocate(top_h2osoi_liq_improad(num_urbanl))
    allocate(top_h2osoi_ice_improad(num_urbanl))

  end subroutine urbanxx_soilFluxes_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilFluxes(num_urbanl, filter_urbanl, &
       num_nolakec, filter_nolakec, &
       num_nolakep, filter_nolakep)
    !
    ! !DESCRIPTION:
    ! Compute soil fluxes (EflxSoilGrnd, QflxEvapGrnd, QflxSubSnow,
    ! QflxDewSnow, QflxDewGrnd) for all five urban surfaces via URBANxx.
    ! Results are written back into the corresponding ELM veg_ef / veg_wf
    ! fields for urban patches.
    !
    use ColumnType     , only : col_pp
    use ColumnDataType  , only : col_ws
    use column_varcon  , only : icol_roof, icol_road_imperv, icol_road_perv, &
                                icol_sunwall, icol_shadewall
    !
    implicit none
    !
    integer(c_int) , intent(in) :: num_urbanl
    integer        , intent(in) :: filter_urbanl(:)   ! urban landunit filter
    integer(c_int) , intent(in) :: num_nolakec
    integer        , intent(in) :: filter_nolakec(:)  ! no-lake column filter
    integer(c_int) , intent(in) :: num_nolakep
    integer        , intent(in) :: filter_nolakep(:)  ! no-lake patch filter
    !
    integer(c_int)  :: status
    integer :: fp, p, c, fc
    integer :: idx_roof, idx_improad, idx_perroad, idx_sunwall, idx_shadwall
    integer :: j  ! top active layer index (snl(c)+1)

    ! --- Pack top-layer soil water for roof and impervious road ---
    ! h2osoi_liq(c, snl(c)+1) and h2osoi_ice(c, snl(c)+1) correspond to
    ! TopH2OSoiLiq and TopH2OSoiIce in URBANxx for roof and impervious road.
    idx_roof    = 0
    idx_improad = 0
    do fc = 1, num_nolakec
      c = filter_nolakec(fc)
      j = col_pp%snl(c) + 1
      if (col_pp%itype(c) == icol_roof) then
        idx_roof = idx_roof + 1
        top_h2osoi_liq_roof(idx_roof) = col_ws%h2osoi_liq(c, j)
        top_h2osoi_ice_roof(idx_roof) = col_ws%h2osoi_ice(c, j)
      else if (col_pp%itype(c) == icol_road_imperv) then
        idx_improad = idx_improad + 1
        top_h2osoi_liq_improad(idx_improad) = col_ws%h2osoi_liq(c, j)
        top_h2osoi_ice_improad(idx_improad) = col_ws%h2osoi_ice(c, j)
      end if
    end do

    ! --- Send top-layer soil water to URBANxx ---
    call UrbanSetTopH2OSoiLiqRoof(urbanxx, c_loc(top_h2osoi_liq_roof), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetTopH2OSoiIceRoof(urbanxx, c_loc(top_h2osoi_ice_roof), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetTopH2OSoiLiqImperviousRoad(urbanxx, c_loc(top_h2osoi_liq_improad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetTopH2OSoiIceImperviousRoad(urbanxx, c_loc(top_h2osoi_ice_improad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! --- Run URBANxx soil fluxes ---
    call UrbanComputeSoilFluxes(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! --- Extract EflxSoilGrnd for all five surfaces ---
    call UrbanGetEflxSoilGrndRoof(urbanxx, c_loc(eflx_soil_grnd_roof), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetEflxSoilGrndImperviousRoad(urbanxx, c_loc(eflx_soil_grnd_improad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetEflxSoilGrndPerviousRoad(urbanxx, c_loc(eflx_soil_grnd_perroad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetEflxSoilGrndSunlitWall(urbanxx, c_loc(eflx_soil_grnd_sunwall), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetEflxSoilGrndShadedWall(urbanxx, c_loc(eflx_soil_grnd_shadwall), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! --- Extract QflxEvapGrnd (roof, imperv, pervious road) ---
    call UrbanGetQflxEvapGrndRoof(urbanxx, c_loc(qflx_evap_grnd_roof), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxEvapGrndImperviousRoad(urbanxx, c_loc(qflx_evap_grnd_improad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxEvapGrndPerviousRoad(urbanxx, c_loc(qflx_evap_grnd_perroad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! --- Extract QflxSubSnow (roof, imperv, pervious road) ---
    call UrbanGetQflxSubSnowRoof(urbanxx, c_loc(qflx_sub_snow_roof), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxSubSnowImperviousRoad(urbanxx, c_loc(qflx_sub_snow_improad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxSubSnowPerviousRoad(urbanxx, c_loc(qflx_sub_snow_perroad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! --- Extract QflxDewSnow (roof, imperv, pervious road) ---
    call UrbanGetQflxDewSnowRoof(urbanxx, c_loc(qflx_dew_snow_roof), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxDewSnowImperviousRoad(urbanxx, c_loc(qflx_dew_snow_improad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxDewSnowPerviousRoad(urbanxx, c_loc(qflx_dew_snow_perroad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! --- Extract QflxDewGrnd (roof, imperv, pervious road) ---
    call UrbanGetQflxDewGrndRoof(urbanxx, c_loc(qflx_dew_grnd_roof), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxDewGrndImperviousRoad(urbanxx, c_loc(qflx_dew_grnd_improad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxDewGrndPerviousRoad(urbanxx, c_loc(qflx_dew_grnd_perroad), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine urbanxx_soilFluxes

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilFluxes_check(num_urbanl, filter_urbanl, &
       num_nolakep, filter_nolakep)
    !
    ! !DESCRIPTION:
    ! Compare URBANxx soil flux results with ELM values.
    ! Calls endrun on mismatch exceeding tolerance.
    !
    use ColumnType    , only : col_pp
    use column_varcon , only : icol_roof, icol_road_imperv, icol_road_perv, &
                               icol_sunwall, icol_shadewall
    !
    implicit none
    !
    integer(c_int) , intent(in) :: num_urbanl
    integer        , intent(in) :: filter_urbanl(:)   ! urban landunit filter
    integer(c_int) , intent(in) :: num_nolakep
    integer        , intent(in) :: filter_nolakep(:)  ! no-lake patch filter
    !
    integer :: fp, p, c
    integer :: idx_roof, idx_improad, idx_perroad, idx_sunwall, idx_shadwall
    real(r8) :: max_error, max_rel_error
    real(r8) :: err_eflx, err_qflx_evap, err_qflx_sub, err_qflx_dsn, err_qflx_dg
    real(r8) :: urb_eflx, urb_qflx_evap, urb_qflx_sub, urb_qflx_dsn, urb_qflx_dg

    max_error     = 0._r8
    max_rel_error = 0._r8
    idx_roof     = 0
    idx_improad  = 0
    idx_perroad  = 0
    idx_sunwall  = 0
    idx_shadwall = 0

    do fp = 1, num_nolakep
      p = filter_nolakep(fp)
      c = veg_pp%column(p)

      select case (col_pp%itype(c))
      case (icol_roof)
        idx_roof = idx_roof + 1
        urb_eflx      = eflx_soil_grnd_roof(idx_roof)
        urb_qflx_evap = qflx_evap_grnd_roof(idx_roof)
        urb_qflx_sub  = qflx_sub_snow_roof(idx_roof)
        urb_qflx_dsn  = qflx_dew_snow_roof(idx_roof)
        urb_qflx_dg   = qflx_dew_grnd_roof(idx_roof)

      case (icol_road_imperv)
        idx_improad = idx_improad + 1
        urb_eflx      = eflx_soil_grnd_improad(idx_improad)
        urb_qflx_evap = qflx_evap_grnd_improad(idx_improad)
        urb_qflx_sub  = qflx_sub_snow_improad(idx_improad)
        urb_qflx_dsn  = qflx_dew_snow_improad(idx_improad)
        urb_qflx_dg   = qflx_dew_grnd_improad(idx_improad)

      case (icol_road_perv)
        idx_perroad = idx_perroad + 1
        urb_eflx      = eflx_soil_grnd_perroad(idx_perroad)
        urb_qflx_evap = qflx_evap_grnd_perroad(idx_perroad)
        urb_qflx_sub  = qflx_sub_snow_perroad(idx_perroad)
        urb_qflx_dsn  = qflx_dew_snow_perroad(idx_perroad)
        urb_qflx_dg   = qflx_dew_grnd_perroad(idx_perroad)

      case (icol_sunwall)
        idx_sunwall = idx_sunwall + 1
        urb_eflx      = eflx_soil_grnd_sunwall(idx_sunwall)
        urb_qflx_evap = 0._r8
        urb_qflx_sub  = 0._r8
        urb_qflx_dsn  = 0._r8
        urb_qflx_dg   = 0._r8

      case (icol_shadewall)
        idx_shadwall = idx_shadwall + 1
        urb_eflx      = eflx_soil_grnd_shadwall(idx_shadwall)
        urb_qflx_evap = 0._r8
        urb_qflx_sub  = 0._r8
        urb_qflx_dsn  = 0._r8
        urb_qflx_dg   = 0._r8

      case default
        cycle

      end select

      err_eflx      = abs(veg_ef%eflx_soil_grnd(p) - urb_eflx)
      err_qflx_evap = abs(veg_wf%qflx_evap_grnd(p) - urb_qflx_evap)
      err_qflx_sub  = abs(veg_wf%qflx_sub_snow(p)  - urb_qflx_sub)
      err_qflx_dsn  = abs(veg_wf%qflx_dew_snow(p)  - urb_qflx_dsn)
      err_qflx_dg   = abs(veg_wf%qflx_dew_grnd(p)  - urb_qflx_dg)

      max_error = max(max_error, err_eflx, err_qflx_evap, err_qflx_sub, &
                      err_qflx_dsn, err_qflx_dg)
      max_rel_error = max(max_rel_error, &
           err_eflx      / max(abs(veg_ef%eflx_soil_grnd(p)), 1.0e-20_r8), &
           err_qflx_evap / max(abs(veg_wf%qflx_evap_grnd(p)), 1.0e-20_r8), &
           err_qflx_sub  / max(abs(veg_wf%qflx_sub_snow(p)),  1.0e-20_r8), &
           err_qflx_dsn  / max(abs(veg_wf%qflx_dew_snow(p)),  1.0e-20_r8), &
           err_qflx_dg   / max(abs(veg_wf%qflx_dew_grnd(p)),  1.0e-20_r8))

      if (err_eflx > 1.0e-10_r8 .or. err_qflx_evap > 1.0e-10_r8 .or. &
          err_qflx_sub > 1.0e-10_r8 .or. err_qflx_dsn > 1.0e-10_r8 .or. &
          err_qflx_dg  > 1.0e-10_r8) then
        write(iulog,*) 'ERROR: Soil flux mismatch between ELM and URBANxx at p=', p, ' c=', c
        write(iulog,*) '  col_type=', col_pp%itype(c)
        write(iulog,*) '  err_eflx_soil_grnd = ', err_eflx
        write(iulog,*) '  err_qflx_evap_grnd = ', err_qflx_evap
        write(iulog,*) '  err_qflx_sub_snow  = ', err_qflx_sub
        write(iulog,*) '  err_qflx_dew_snow  = ', err_qflx_dsn
        write(iulog,*) '  err_qflx_dew_grnd  = ', err_qflx_dg
        write(iulog,*) '  ELM    eflx_soil_grnd=', veg_ef%eflx_soil_grnd(p)
        write(iulog,*) '  URBANxx eflx_soil_grnd=', urb_eflx
        call exit(0)
      end if

    end do

    write(iulog,*) 'Max error in soil ground heat flux  : ', max_error, ' (rel: ', max_rel_error, ')'

  end subroutine urbanxx_soilFluxes_check

end module UrbanxxSoilFluxesMod
