module UrbanxxWaterTableMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Set water table boundary conditions and compute WaterTable physics
  ! for pervious road columns in the Urban++ model.
  !
  ! Ports ELM's WaterTable subroutine (SoilHydrologyMod.F90:746-1082)
  ! for the pervious road column only.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use urban_kokkos_interface
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use ColumnDataType       , only : col_ws, col_wf
  use UrbanxxInstanceMod   , only : urbanxx

  implicit none

  private

  ! Persistent input buffers (allocated once in init)
  real(c_double), allocatable, target :: wa_in(:)
  real(c_double), allocatable, target :: frac_h2osfc_in(:)
  real(c_double), allocatable, target :: dew_grnd_in(:)
  real(c_double), allocatable, target :: dew_snow_in(:)
  real(c_double), allocatable, target :: sub_snow_in(:)
  real(c_double), allocatable, target :: qcharge_in(:)

  ! Persistent output buffers (allocated once in init)
  real(c_double), allocatable, target, public :: out_zwt(:)
  real(c_double), allocatable, target, public :: out_wa(:)
  real(c_double), allocatable, target, public :: out_zwt_perched(:)
  real(c_double), allocatable, target, public :: out_qflx_sub_snow(:)
  real(c_double), allocatable, target, public :: out_qflx_drain(:)
  real(c_double), allocatable, target, public :: out_qflx_rsub_sat(:)
  ! 2D flattened: (num_urbanl * nlevgrnd)
  real(c_double), allocatable, target, public :: out_h2osoi_liq(:)
  real(c_double), allocatable, target, public :: out_h2osoi_ice(:)

  ! -----------------------------------------------------------------------
  ! Dew condensation buffers (roof and impervious road)
  ! -----------------------------------------------------------------------
  ! Input buffers
  real(c_double), allocatable, target :: dew_grnd_roof_in(:)
  real(c_double), allocatable, target :: dew_snow_roof_in(:)
  real(c_double), allocatable, target :: sub_snow_roof_in(:)
  real(c_double), allocatable, target :: top_liq_roof_in(:)
  real(c_double), allocatable, target :: top_ice_roof_in(:)

  real(c_double), allocatable, target :: dew_grnd_imperv_in(:)
  real(c_double), allocatable, target :: dew_snow_imperv_in(:)
  real(c_double), allocatable, target :: sub_snow_imperv_in(:)
  real(c_double), allocatable, target :: top_liq_imperv_in(:)
  real(c_double), allocatable, target :: top_ice_imperv_in(:)

  ! Output buffers
  real(c_double), allocatable, target, public :: out_top_liq_roof(:)
  real(c_double), allocatable, target, public :: out_top_ice_roof(:)
  real(c_double), allocatable, target, public :: out_qflx_sub_snow_roof(:)

  real(c_double), allocatable, target, public :: out_top_liq_imperv(:)
  real(c_double), allocatable, target, public :: out_top_ice_imperv(:)
  real(c_double), allocatable, target, public :: out_qflx_sub_snow_imperv(:)

  public :: urbanxx_waterTable_init
  public :: urbanxx_waterTable
  public :: urbanxx_waterTable_check
  public :: urbanxx_dewCondensation
  public :: urbanxx_dewCondensation_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_waterTable_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for water table computation.
    ! Called once during initialization.
    !
    use elm_varpar, only : nlevgrnd
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int) :: totalSize

    totalSize = num_urbanl * nlevgrnd

    ! Input buffers — 1D
    allocate(wa_in(num_urbanl))
    allocate(frac_h2osfc_in(num_urbanl))
    allocate(dew_grnd_in(num_urbanl))
    allocate(dew_snow_in(num_urbanl))
    allocate(sub_snow_in(num_urbanl))
    allocate(qcharge_in(num_urbanl))

    ! Output buffers — 1D
    allocate(out_zwt(num_urbanl))
    allocate(out_wa(num_urbanl))
    allocate(out_zwt_perched(num_urbanl))
    allocate(out_qflx_sub_snow(num_urbanl))
    allocate(out_qflx_drain(num_urbanl))
    allocate(out_qflx_rsub_sat(num_urbanl))

    ! Output buffers — 2D flattened
    allocate(out_h2osoi_liq(totalSize))
    allocate(out_h2osoi_ice(totalSize))

    ! Dew condensation input buffers (1D)
    allocate(dew_grnd_roof_in(num_urbanl))
    allocate(dew_snow_roof_in(num_urbanl))
    allocate(sub_snow_roof_in(num_urbanl))
    allocate(top_liq_roof_in(num_urbanl))
    allocate(top_ice_roof_in(num_urbanl))

    allocate(dew_grnd_imperv_in(num_urbanl))
    allocate(dew_snow_imperv_in(num_urbanl))
    allocate(sub_snow_imperv_in(num_urbanl))
    allocate(top_liq_imperv_in(num_urbanl))
    allocate(top_ice_imperv_in(num_urbanl))

    ! Dew condensation output buffers (1D)
    allocate(out_top_liq_roof(num_urbanl))
    allocate(out_top_ice_roof(num_urbanl))
    allocate(out_qflx_sub_snow_roof(num_urbanl))

    allocate(out_top_liq_imperv(num_urbanl))
    allocate(out_top_ice_imperv(num_urbanl))
    allocate(out_qflx_sub_snow_imperv(num_urbanl))

  end subroutine urbanxx_waterTable_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_waterTable(num_urbanl, num_urbanc, filter_urbanc, &
                                soilhydrology_vars, dtime)
    !
    ! !DESCRIPTION:
    ! Set water table inputs for pervious road columns, call
    ! UrbanComputeWaterTable, and unpack outputs into ELM arrays.
    !
    use ColumnType           , only : col_pp
    use column_varcon        , only : icol_road_perv
    use elm_varpar           , only : nlevgrnd
    use SoilHydrologyType    , only : soilhydrology_type
    use abortutils           , only : endrun
    use shr_log_mod          , only : errmsg => shr_log_errmsg
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer       , intent(in) :: filter_urbanc(:)
    type(soilhydrology_type), intent(inout) :: soilhydrology_vars
    real(r8)      , intent(in) :: dtime
    !
    ! !LOCAL VARIABLES:
    integer(c_int)               :: status
    integer                      :: fc, c, j, idx, idx_perv, nlevbed
    integer(c_int)               :: totalSize
    integer(c_int), dimension(2) :: size2D
    logical(c_bool)              :: isLayoutLeft

    associate(                                                               &
         h2osoi_liq     => col_ws%h2osoi_liq                              , & ! Output: liquid water [kg/m2]
         h2osoi_ice     => col_ws%h2osoi_ice                              , & ! Output: ice lens [kg/m2]
         frac_h2osfc    => col_ws%frac_h2osfc                             , & ! Input:  fraction surface covered by ponded water [-]
         qflx_dew_grnd  => col_wf%qflx_dew_grnd                          , & ! Input:  ground dew flux [mm H2O/s]
         qflx_dew_snow  => col_wf%qflx_dew_snow                          , & ! Input:  dew added to snow [mm H2O/s]
         qflx_sub_snow  => col_wf%qflx_sub_snow                          , & ! In/Out: sublimation from ice [mm H2O/s]
         qflx_drain     => col_wf%qflx_drain                             , & ! Output: sub-surface drainage [mm H2O/s]
         qflx_rsub_sat  => col_wf%qflx_rsub_sat                          , & ! Output: saturation excess runoff [mm H2O/s]
         zwt_col        => soilhydrology_vars%zwt_col                     , & ! In/Out: water table depth [m]
         wa_col         => soilhydrology_vars%wa_col                      , & ! In/Out: aquifer water [mm]
         zwt_perched_col => soilhydrology_vars%zwt_perched_col            , & ! Output: perched water table [m]
         qcharge_col    => soilhydrology_vars%qcharge_col                 , & ! Input:  aquifer recharge rate [mm/s]
         nlev2bed       => col_pp%nlevbed                                   & ! Input:  number of layers to bedrock
         )

      ! --------------------------------------------------------
      ! Pack 1D input buffers for pervious road columns
      ! --------------------------------------------------------
      idx_perv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          wa_in(idx_perv)         = wa_col(c)
          frac_h2osfc_in(idx_perv) = frac_h2osfc(c)
          dew_grnd_in(idx_perv)   = qflx_dew_grnd(c)
          dew_snow_in(idx_perv)   = qflx_dew_snow(c)
          sub_snow_in(idx_perv)   = qflx_sub_snow(c)
          qcharge_in(idx_perv)    = qcharge_col(c)
        end if
      end do

      ! Set 1D inputs
      call UrbanSetAquiferWaterForPerviousRoad(urbanxx, c_loc(wa_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetFracH2osfcForPerviousRoad(urbanxx, c_loc(frac_h2osfc_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetDewGrndFluxForPerviousRoad(urbanxx, c_loc(dew_grnd_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetDewSnowFluxForPerviousRoad(urbanxx, c_loc(dew_snow_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetSubSnowFluxForPerviousRoad(urbanxx, c_loc(sub_snow_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetQchargeForPerviousRoad(urbanxx, c_loc(qcharge_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! --------------------------------------------------------
      ! Set h2osoi_liq and h2osoi_ice from current ELM values
      ! (UrbanComputeWaterTable is called BEFORE ELM's WaterTable modifies them)
      ! Reuse out_h2osoi_liq/ice as pack buffers; they are overwritten by getters later.
      ! --------------------------------------------------------
      totalSize = num_urbanl * nlevgrnd
      size2D(1) = num_urbanl
      size2D(2) = nlevgrnd
      isLayoutLeft = UrbanKokkosIsLayoutLeft()

      if (isLayoutLeft) then
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
                out_h2osoi_liq(idx) = h2osoi_liq(c, j)
                out_h2osoi_ice(idx) = h2osoi_ice(c, j)
              else
                out_h2osoi_liq(idx) = 0.0_r8
                out_h2osoi_ice(idx) = 0.0_r8
              end if
            end if
          end do
        end do
      else
        idx = 0
        do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col_pp%itype(c) == icol_road_perv) then
            nlevbed = nlev2bed(c)
            do j = 1, nlevgrnd
              idx = idx + 1
              if (j <= nlevbed) then
                out_h2osoi_liq(idx) = h2osoi_liq(c, j)
                out_h2osoi_ice(idx) = h2osoi_ice(c, j)
              else
                out_h2osoi_liq(idx) = 0.0_r8
                out_h2osoi_ice(idx) = 0.0_r8
              end if
            end do
          end if
        end do
      end if

      call UrbanSetSoilLiquidWaterForPerviousRoad(urbanxx, c_loc(out_h2osoi_liq), &
           size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetSoilIceContentForPerviousRoad(urbanxx, c_loc(out_h2osoi_ice), &
           size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! --------------------------------------------------------
      ! Compute
      ! --------------------------------------------------------
      call UrbanComputeWaterTable(urbanxx, dtime, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! --------------------------------------------------------
      ! Retrieve 1D outputs
      ! --------------------------------------------------------
      call UrbanGetWaterTableDepthPerviousRoad(urbanxx, c_loc(out_zwt), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetAquiferWaterPerviousRoad(urbanxx, c_loc(out_wa), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetZwtPerchedPerviousRoad(urbanxx, c_loc(out_zwt_perched), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetSubSnowFluxPerviousRoad(urbanxx, c_loc(out_qflx_sub_snow), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetDrainFluxPerviousRoad(urbanxx, c_loc(out_qflx_drain), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetRsubSatPerviousRoad(urbanxx, c_loc(out_qflx_rsub_sat), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! --------------------------------------------------------
      ! Retrieve 2D outputs (h2osoi_liq, h2osoi_ice)
      ! --------------------------------------------------------
      isLayoutLeft = UrbanKokkosIsLayoutLeft()

      call UrbanGetSoilLiquidWaterPerviousRoad(urbanxx, c_loc(out_h2osoi_liq), &
           size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetSoilIceContentPerviousRoad(urbanxx, c_loc(out_h2osoi_ice), &
           size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    end associate

  end subroutine urbanxx_waterTable

  !-----------------------------------------------------------------------
  subroutine urbanxx_waterTable_check(num_urbanl, num_urbanc, filter_urbanc, &
                                        soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Compare UrbanXX WaterTable outputs against ELM values for pervious road.
    ! Reports max absolute and relative errors to the log.
    ! Calls endrun if any error exceeds tolerance.
    !
    use ColumnType            , only : col_pp
    use column_varcon         , only : icol_road_perv
    use elm_varpar            , only : nlevgrnd
    use abortutils            , only : endrun
    use shr_log_mod           , only : errmsg => shr_log_errmsg
    use SoilHydrologyType     , only : soilhydrology_type
    use urban_kokkos_interface, only : UrbanKokkosIsLayoutLeft
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int)         , intent(in) :: num_urbanl
    integer(c_int)         , intent(in) :: num_urbanc
    integer                , intent(in) :: filter_urbanc(:)
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    !
    ! !LOCAL VARIABLES:
    ! (outputs are in module-level out_* buffers from the last call to urbanxx_waterTable)
    integer  :: fc, c, j, l, count, idx_perv
    real(r8) :: max_abs_err, max_rel_err, abs_err

    ! Tolerance for endrun
    real(r8), parameter :: tol = 1.0e-6_r8

    associate( &
         h2osoi_liq      => col_ws%h2osoi_liq                       , &
         h2osoi_ice      => col_ws%h2osoi_ice                       , &
         zwt_col         => soilhydrology_vars%zwt_col              , &
         wa_col          => soilhydrology_vars%wa_col               , &
         zwt_perched_col => soilhydrology_vars%zwt_perched_col        &
         )

      ! --- zwt_col ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8
      idx_perv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          abs_err = abs(zwt_col(c) - out_zwt(idx_perv))
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, &
               abs_err / max(abs(zwt_col(c)), 1.0e-20_r8))
        end if
      end do
      write(iulog,*) 'Max error in zwt_col                :', max_abs_err, '  (rel:', max_rel_err, ')'
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_waterTable_check: zwt_col error too large'//errmsg(__FILE__,__LINE__))

      ! --- wa_col ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8
      idx_perv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          abs_err = abs(wa_col(c) - out_wa(idx_perv))
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, &
               abs_err / max(abs(wa_col(c)), 1.0e-20_r8))
        end if
      end do
      write(iulog,*) 'Max error in wa_col                 :', max_abs_err, '  (rel:', max_rel_err, ')'
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_waterTable_check: wa_col error too large'//errmsg(__FILE__,__LINE__))

      ! --- zwt_perched_col ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8
      idx_perv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          abs_err = abs(zwt_perched_col(c) - out_zwt_perched(idx_perv))
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, &
               abs_err / max(abs(zwt_perched_col(c)), 1.0e-20_r8))
        end if
      end do
      write(iulog,*) 'Max error in zwt_perched_col        :', max_abs_err, '  (rel:', max_rel_err, ')'
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_waterTable_check: zwt_perched_col error too large'//errmsg(__FILE__,__LINE__))

      ! --- h2osoi_liq ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8

      if (UrbanKokkosIsLayoutLeft()) then
        count = 0
        do j = 1, nlevgrnd
          idx_perv = 0
          do fc = 1, num_urbanc
            c = filter_urbanc(fc)
            if (col_pp%itype(c) == icol_road_perv) then
              idx_perv = idx_perv + 1
              count = count + 1
              abs_err = abs(h2osoi_liq(c,j) - out_h2osoi_liq(count))
              max_abs_err = max(max_abs_err, abs_err)
              max_rel_err = max(max_rel_err, &
                   abs_err / max(abs(h2osoi_liq(c,j)), 1.0e-20_r8))
            end if
          end do
        end do
      else
        count = 0
        do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col_pp%itype(c) == icol_road_perv) then
            do j = 1, nlevgrnd
              count = count + 1
              abs_err = abs(h2osoi_liq(c,j) - out_h2osoi_liq(count))
              max_abs_err = max(max_abs_err, abs_err)
              max_rel_err = max(max_rel_err, &
                   abs_err / max(abs(h2osoi_liq(c,j)), 1.0e-20_r8))
            end do
          end if
        end do
      end if
      write(iulog,*) 'Max error in h2osoi_liq             :', max_abs_err, '  (rel:', max_rel_err, ')'
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_waterTable_check: h2osoi_liq error too large'//errmsg(__FILE__,__LINE__))

      ! --- h2osoi_ice ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8

      if (UrbanKokkosIsLayoutLeft()) then
        count = 0
        do j = 1, nlevgrnd
          idx_perv = 0
          do fc = 1, num_urbanc
            c = filter_urbanc(fc)
            if (col_pp%itype(c) == icol_road_perv) then
              idx_perv = idx_perv + 1
              count = count + 1
              abs_err = abs(h2osoi_ice(c,j) - out_h2osoi_ice(count))
              max_abs_err = max(max_abs_err, abs_err)
              max_rel_err = max(max_rel_err, &
                   abs_err / max(abs(h2osoi_ice(c,j)), 1.0e-20_r8))
            end if
          end do
        end do
      else
        count = 0
        do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col_pp%itype(c) == icol_road_perv) then
            do j = 1, nlevgrnd
              count = count + 1
              abs_err = abs(h2osoi_ice(c,j) - out_h2osoi_ice(count))
              max_abs_err = max(max_abs_err, abs_err)
              max_rel_err = max(max_rel_err, &
                   abs_err / max(abs(h2osoi_ice(c,j)), 1.0e-20_r8))
            end do
          end if
        end do
      end if
      write(iulog,*) 'Max error in h2osoi_ice             :', max_abs_err, '  (rel:', max_rel_err, ')'
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_waterTable_check: h2osoi_ice error too large'//errmsg(__FILE__,__LINE__))

    end associate

  end subroutine urbanxx_waterTable_check

  !-----------------------------------------------------------------------
  subroutine urbanxx_dewCondensation(num_urbanl, num_urbanc, filter_urbanc, &
                                     soilhydrology_vars, dtime)
    !
    ! !DESCRIPTION:
    ! Set dew/sublimation inputs for roof and impervious road columns,
    ! call UrbanComputeDewCondensationRoofImperviousRoad, and retrieve outputs.
    ! ELM reference: SoilHydrologyMod.F90:1062-1078
    !
    use ColumnType           , only : col_pp
    use column_varcon        , only : icol_roof, icol_road_imperv
    use SoilHydrologyType    , only : soilhydrology_type
    use abortutils           , only : endrun
    use shr_log_mod          , only : errmsg => shr_log_errmsg
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer       , intent(in) :: filter_urbanc(:)
    type(soilhydrology_type), intent(inout) :: soilhydrology_vars
    real(r8)      , intent(in) :: dtime
    !
    ! !LOCAL VARIABLES:
    integer(c_int) :: status
    integer        :: fc, c, l_roof, l_imperv

    associate(                                                            &
         h2osoi_liq    => col_ws%h2osoi_liq                           , & ! In:  liquid [kg/m2]
         h2osoi_ice    => col_ws%h2osoi_ice                           , & ! In:  ice    [kg/m2]
         qflx_dew_grnd => col_wf%qflx_dew_grnd                       , & ! In:  ground dew flux [mm/s]
         qflx_dew_snow => col_wf%qflx_dew_snow                       , & ! In:  dew to snow [mm/s]
         qflx_sub_snow => col_wf%qflx_sub_snow                         & ! In:  sublimation from ice [mm/s]
         )

      ! --------------------------------------------------------
      ! Pack 1D input buffers for roof and impervious road columns
      ! --------------------------------------------------------
      l_roof   = 0
      l_imperv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_roof) then
          l_roof = l_roof + 1
          top_liq_roof_in(l_roof)   = h2osoi_liq(c,1)
          top_ice_roof_in(l_roof)   = h2osoi_ice(c,1)
          dew_grnd_roof_in(l_roof)  = qflx_dew_grnd(c)
          dew_snow_roof_in(l_roof)  = qflx_dew_snow(c)
          sub_snow_roof_in(l_roof)  = qflx_sub_snow(c)
        else if (col_pp%itype(c) == icol_road_imperv) then
          l_imperv = l_imperv + 1
          top_liq_imperv_in(l_imperv)  = h2osoi_liq(c,1)
          top_ice_imperv_in(l_imperv)  = h2osoi_ice(c,1)
          dew_grnd_imperv_in(l_imperv) = qflx_dew_grnd(c)
          dew_snow_imperv_in(l_imperv) = qflx_dew_snow(c)
          sub_snow_imperv_in(l_imperv) = qflx_sub_snow(c)
        end if
      end do

      ! --------------------------------------------------------
      ! Set roof inputs
      ! --------------------------------------------------------
      call UrbanSetTopH2OSoiLiqRoof(urbanxx, c_loc(top_liq_roof_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetTopH2OSoiIceRoof(urbanxx, c_loc(top_ice_roof_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetQflxDewGrndRoof(urbanxx, c_loc(dew_grnd_roof_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetQflxDewSnowRoof(urbanxx, c_loc(dew_snow_roof_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetQflxSubSnowRoof(urbanxx, c_loc(sub_snow_roof_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! --------------------------------------------------------
      ! Set impervious road inputs
      ! --------------------------------------------------------
      call UrbanSetTopH2OSoiLiqImperviousRoad(urbanxx, c_loc(top_liq_imperv_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetTopH2OSoiIceImperviousRoad(urbanxx, c_loc(top_ice_imperv_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetQflxDewGrndImperviousRoad(urbanxx, c_loc(dew_grnd_imperv_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetQflxDewSnowImperviousRoad(urbanxx, c_loc(dew_snow_imperv_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetQflxSubSnowImperviousRoad(urbanxx, c_loc(sub_snow_imperv_in), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! --------------------------------------------------------
      ! Compute
      ! --------------------------------------------------------
      call UrbanComputeDewCondensationRoofImperviousRoad(urbanxx, dtime, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! --------------------------------------------------------
      ! Retrieve outputs
      ! --------------------------------------------------------
      call UrbanGetTopH2OSoiLiqRoof(urbanxx, c_loc(out_top_liq_roof), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetTopH2OSoiIceRoof(urbanxx, c_loc(out_top_ice_roof), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetQflxSubSnowRoof(urbanxx, c_loc(out_qflx_sub_snow_roof), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetTopH2OSoiLiqImperviousRoad(urbanxx, c_loc(out_top_liq_imperv), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetTopH2OSoiIceImperviousRoad(urbanxx, c_loc(out_top_ice_imperv), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetQflxSubSnowImperviousRoad(urbanxx, c_loc(out_qflx_sub_snow_imperv), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    end associate

    ! No output values are written back into ELM data structures.

  end subroutine urbanxx_dewCondensation

  !-----------------------------------------------------------------------
  subroutine urbanxx_dewCondensation_check(num_urbanl, num_urbanc, &
                                           filter_urbanc, soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! Compare URBANxx dew condensation outputs against ELM values for
    ! roof and impervious road columns. Calls endrun on any mismatch.
    !
    use ColumnType   , only : col_pp
    use column_varcon, only : icol_roof, icol_road_imperv
    use abortutils   , only : endrun
    use shr_log_mod  , only : errmsg => shr_log_errmsg
    use SoilHydrologyType , only : soilhydrology_type
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int)          , intent(in) :: num_urbanl
    integer(c_int)          , intent(in) :: num_urbanc
    integer                 , intent(in) :: filter_urbanc(:)
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, l_roof, l_imperv
    real(r8) :: max_abs_err, max_rel_err, abs_err
    real(r8), parameter :: tol = 1.0e-6_r8

    associate( &
         h2osoi_liq => col_ws%h2osoi_liq , &
         h2osoi_ice => col_ws%h2osoi_ice   &
         )

      ! --- h2osoi_liq (roof) ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8
      l_roof = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_roof) then
          l_roof = l_roof + 1
          abs_err = abs(h2osoi_liq(c,1) - out_top_liq_roof(l_roof))
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, &
               abs_err / max(abs(h2osoi_liq(c,1)), 1.0e-20_r8))
        end if
      end do
      write(iulog,*) 'Max error in h2osoi_liq (roof)      :', max_abs_err, '  (rel:', max_rel_err, ')'
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_dewCondensation_check: h2osoi_liq roof error too large'//errmsg(__FILE__,__LINE__))

      ! --- h2osoi_ice (roof) ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8
      l_roof = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_roof) then
          l_roof = l_roof + 1
          abs_err = abs(h2osoi_ice(c,1) - out_top_ice_roof(l_roof))
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, &
               abs_err / max(abs(h2osoi_ice(c,1)), 1.0e-20_r8))
        end if
      end do
      write(iulog,*) 'Max error in h2osoi_ice (roof)      :', max_abs_err, '  (rel:', max_rel_err, ')'
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_dewCondensation_check: h2osoi_ice roof error too large'//errmsg(__FILE__,__LINE__))

      ! --- h2osoi_liq (impervious road) ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8
      l_imperv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_imperv) then
          l_imperv = l_imperv + 1
          abs_err = abs(h2osoi_liq(c,1) - out_top_liq_imperv(l_imperv))
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, &
               abs_err / max(abs(h2osoi_liq(c,1)), 1.0e-20_r8))
        end if
      end do
      write(iulog,*) 'Max error in h2osoi_liq (imperv)    :', max_abs_err, '  (rel:', max_rel_err, ')'
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_dewCondensation_check: h2osoi_liq imperv road error too large'//errmsg(__FILE__,__LINE__))

      ! --- h2osoi_ice (impervious road) ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8
      l_imperv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_imperv) then
          l_imperv = l_imperv + 1
          abs_err = abs(h2osoi_ice(c,1) - out_top_ice_imperv(l_imperv))
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, &
               abs_err / max(abs(h2osoi_ice(c,1)), 1.0e-20_r8))
        end if
      end do
      write(iulog,*) 'Max error in h2osoi_ice (imperv)    :', max_abs_err, '  (rel:', max_rel_err, ')'
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_dewCondensation_check: h2osoi_ice imperv road error too large'//errmsg(__FILE__,__LINE__))

    end associate

  end subroutine urbanxx_dewCondensation_check

end module UrbanxxWaterTableMod
