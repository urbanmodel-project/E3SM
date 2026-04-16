module UrbanxxDrainageMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Set drainage inputs and compute Drainage physics for pervious road
  ! columns in the Urban++ model.
  !
  ! Ports ELM's Drainage subroutine (SoilHydrologyMod.F90:1085-1773)
  ! for the pervious road column only.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use urban_kokkos_interface
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use ColumnDataType       , only : col_wf
  use UrbanxxInstanceMod   , only : urbanxx

  implicit none

  private

  ! Number of urban landunits (set at initialization, used by drainage routines)
  integer(c_int), save :: module_num_urbanl = 0

  ! Persistent input buffers (allocated once in init)
  real(c_double), allocatable, target :: hkdepth_in(:)
  real(c_double), allocatable, target :: topo_slope_in(:)

  ! Persistent output buffers (allocated once in init)
  real(c_double), allocatable, target, public :: out_qflx_drain(:)
  real(c_double), allocatable, target, public :: out_qflx_rsub_sat(:)

  public :: urbanxx_drainage_init
  public :: urbanxx_drainage
  public :: urbanxx_drainage_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_drainage_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers and send init-time constants to URBANxx.
    ! Called once during initialization.
    !
    use pftvarcon        , only : rsub_top_globalmax
    use elm_varcon       , only : pondmx, watmin, e_ice
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    !
    integer(c_int) :: status

    module_num_urbanl = num_urbanl

    ! Allocate input buffers
    allocate(hkdepth_in(num_urbanl))
    allocate(topo_slope_in(num_urbanl))

    ! Allocate output buffers
    allocate(out_qflx_drain(num_urbanl))
    allocate(out_qflx_rsub_sat(num_urbanl))

    ! Send constants to URBANxx
    call UrbanSetRsubTopGlobalMax(urbanxx, real(rsub_top_globalmax, c_double), status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetPondmax(urbanxx, real(pondmx, c_double), status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetWatmin(urbanxx, real(watmin, c_double), status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetEice(urbanxx, real(e_ice, c_double), status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine urbanxx_drainage_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_drainage(num_urbanc, filter_urbanc, &
                               soilhydrology_vars, dtime)
    !
    ! !DESCRIPTION:
    ! Set drainage inputs for pervious road columns, call
    ! UrbanComputeDrainage, and unpack outputs into module buffers.
    !
    use ColumnType           , only : col_pp
    use column_varcon        , only : icol_road_perv
    use SoilHydrologyType    , only : soilhydrology_type
    use abortutils           , only : endrun
    use shr_log_mod          , only : errmsg => shr_log_errmsg
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent(in) :: num_urbanc
    integer, intent(in) :: filter_urbanc(:)
    type(soilhydrology_type), intent(inout) :: soilhydrology_vars
    real(r8)      , intent(in) :: dtime
    !
    ! !LOCAL VARIABLES:
    integer(c_int) :: status
    integer        :: fc, c, idx_perv

    associate( &
         hkdepth    => soilhydrology_vars%hkdepth_col , &
         topo_slope => col_pp%topo_slope                &
         )

      ! --------------------------------------------------------
      ! Pack 1D input buffers for pervious road columns
      ! --------------------------------------------------------
      idx_perv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          hkdepth_in(idx_perv)    = hkdepth(c)
          topo_slope_in(idx_perv) = topo_slope(c)
        end if
      end do

      ! Set inputs
      call UrbanSetHkDepthForPerviousRoad(urbanxx, c_loc(hkdepth_in), &
           module_num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetTopoSlopeForPerviousRoad(urbanxx, c_loc(topo_slope_in), &
           module_num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! --------------------------------------------------------
      ! Compute
      ! --------------------------------------------------------
      call UrbanComputeDrainage(urbanxx, dtime, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! --------------------------------------------------------
      ! Retrieve outputs
      ! --------------------------------------------------------
      call UrbanGetDrainFluxPerviousRoad(urbanxx, c_loc(out_qflx_drain), &
           module_num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetRsubSatPerviousRoad(urbanxx, c_loc(out_qflx_rsub_sat), &
           module_num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    end associate

  end subroutine urbanxx_drainage

  !-----------------------------------------------------------------------
  subroutine urbanxx_drainage_check(num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Compare URBANxx Drainage outputs against ELM values for pervious road.
    ! Reports max absolute and relative errors to the log.
    ! Calls endrun if any error exceeds tolerance.
    !
    use ColumnType        , only : col_pp
    use column_varcon     , only : icol_road_perv
    use abortutils        , only : endrun
    use shr_log_mod       , only : errmsg => shr_log_errmsg
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer, intent(in) :: num_urbanc
    integer, intent(in) :: filter_urbanc(:)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, idx_perv
    real(r8) :: max_abs_err, max_rel_err, abs_err

    ! Tolerance for endrun
    real(r8), parameter :: tol = 1.0e-6_r8

    associate( &
         qflx_drain    => col_wf%qflx_drain    , &
         qflx_rsub_sat => col_wf%qflx_rsub_sat   &
         )

      ! --- qflx_drain ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8
      idx_perv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          abs_err = abs(qflx_drain(c) - out_qflx_drain(idx_perv))
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, &
               abs_err / max(abs(qflx_drain(c)), 1.0e-20_r8))
        end if
      end do
      if (masterproc) then
         write(iulog,*) 'Max error in qflx_drain             :', &
             max_abs_err, '  (rel:', max_rel_err, ')'
      end if
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_drainage_check: qflx_drain error too large'// &
           errmsg(__FILE__,__LINE__))

      ! --- qflx_rsub_sat ---
      max_abs_err = 0._r8
      max_rel_err = 0._r8
      idx_perv = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          abs_err = abs(qflx_rsub_sat(c) - out_qflx_rsub_sat(idx_perv))
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, &
               abs_err / max(abs(qflx_rsub_sat(c)), 1.0e-20_r8))
        end if
      end do
      if (masterproc) then
        write(iulog,*) 'Max error in qflx_rsub_sat           :', &
             max_abs_err, '  (rel:', max_rel_err, ')'
      end if
      if (max_abs_err > tol) call endrun( &
           msg='urbanxx_drainage_check: qflx_rsub_sat error too large'// &
           errmsg(__FILE__,__LINE__))

    end associate

  end subroutine urbanxx_drainage_check

end module UrbanxxDrainageMod
