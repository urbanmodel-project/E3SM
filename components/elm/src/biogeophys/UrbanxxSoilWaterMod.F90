module UrbanxxSoilWaterMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Set soil water boundary conditions and compute hydrology
  ! for urban areas in the Urban++ model.
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

  ! Persistent output buffers (allocated once in init)
  ! 2D: (num_urbanl * nlevgrnd) flattened
  real(c_double) , allocatable, target, public :: out_h2osoi_liq(:)
  real(c_double) , allocatable, target, public :: out_h2osoi_vol(:)
  ! 1D: (num_urbanl)
  real(c_double) , allocatable, target, public :: out_qcharge(:)
  real(c_double) , allocatable, target, public :: out_qflx_deficit(:)

  public :: urbanxx_soilWater_init
  public :: urbanxx_soilWater
  public :: urbanxx_soilWater_check


contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilWater_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for hydrology computation.
    ! Called once during initialization.
    !
    use elm_varpar, only : nlevgrnd
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int) :: totalSize

    totalSize = num_urbanl * nlevgrnd

    ! Output buffers — 2D flattened
    allocate(out_h2osoi_liq(totalSize))
    allocate(out_h2osoi_vol(totalSize))

    ! Output buffers — 1D
    allocate(out_qcharge(num_urbanl))
    allocate(out_qflx_deficit(num_urbanl))

  end subroutine urbanxx_soilWater_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilWater(num_urbanl, dtime)
    !
    ! !DESCRIPTION:
    ! Compute soil hydrology for urban pervious road columns.
    ! Root-fraction-weighted transpiration is computed internally in URBANxx
    ! using Rootr (seeded at initialization from ELM's rootr_col) and
    ! QflxTranEvap (computed each timestep by UrbanComputeSurfaceFluxes).
    !
    use elm_varpar   , only : nlevgrnd
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int), intent(in) :: num_urbanl
    real(r8)      , intent(in) :: dtime                ! time step (s)
    !
    ! !LOCAL VARIABLES:
    integer(c_int)               :: status
    integer(c_int)               :: totalSize
    integer(c_int), dimension(2) :: size2D

    totalSize = num_urbanl * nlevgrnd
    size2D(1) = num_urbanl
    size2D(2) = nlevgrnd

    call UrbanComputeHydrology(urbanxx, dtime, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! Extract hydrology outputs from UrbanXX
      ! 2D outputs: soil liquid water and volumetric water (pervious road)
      call UrbanGetSoilLiquidWaterPerviousRoad(urbanxx, c_loc(out_h2osoi_liq), size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetSoilVolumetricWaterPerviousRoad(urbanxx, c_loc(out_h2osoi_vol), size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      ! 1D outputs: aquifer recharge rate and water deficit flux
      call UrbanGetAquiferRechargeRatePerviousRoad(urbanxx, c_loc(out_qcharge), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)
      call UrbanGetWaterDeficitFluxPerviousRoad(urbanxx, c_loc(out_qflx_deficit), num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine urbanxx_soilWater

  !-----------------------------------------------------------------------
  subroutine urbanxx_soilWater_check(num_urbanl, num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Compare UrbanXX soil water outputs against ELM values for pervious road.
    !
    use ColumnType     , only : col_pp
    use column_varcon  , only : icol_road_perv
    use elm_varpar     , only : nlevgrnd
    use urban_kokkos_interface , only : UrbanKokkosIsLayoutLeft
    use elm_time_manager       , only : get_nstep
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer       , intent(in) :: filter_urbanc(:)
    !
    ! !LOCAL VARIABLES:
    real(c_double), allocatable, target :: h2osoi_liq_2d(:,:)
    real(c_double), allocatable, target :: h2osoi_vol_2d(:,:)
    integer  :: fc, c, j, l, count, idx_perv
    integer  :: max_err_c, max_err_j, max_err_idx_perv
    real(r8) :: max_error_liq, max_error_vol, max_rel_error_liq
    real(r8) :: diff, max_err_elm, max_err_uxx

    associate( &
         h2osoi_liq => col_ws%h2osoi_liq , &
         h2osoi_vol => col_ws%h2osoi_vol   &
         )

      allocate(h2osoi_liq_2d(num_urbanl, nlevgrnd))
      allocate(h2osoi_vol_2d(num_urbanl, nlevgrnd))

      ! Reshape flattened output arrays into 2D based on Kokkos layout
      if (UrbanKokkosIsLayoutLeft()) then
        count = 0
        do j = 1, nlevgrnd
          do l = 1, num_urbanl
            count = count + 1
            h2osoi_liq_2d(l,j) = out_h2osoi_liq(count)
            h2osoi_vol_2d(l,j) = out_h2osoi_vol(count)
          end do
        end do
      else
        count = 0
        do l = 1, num_urbanl
          do j = 1, nlevgrnd
            count = count + 1
            h2osoi_liq_2d(l,j) = out_h2osoi_liq(count)
            h2osoi_vol_2d(l,j) = out_h2osoi_vol(count)
          end do
        end do
      end if

      max_error_liq     = 0._r8
      max_error_vol     = 0._r8
      max_rel_error_liq = 0._r8
      max_err_c         = -1
      max_err_j         = -1
      max_err_idx_perv  = -1
      max_err_elm       = 0._r8
      max_err_uxx       = 0._r8
      idx_perv = 0

      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          do j = 1, nlevgrnd
            diff              = abs(h2osoi_liq(c,j) - h2osoi_liq_2d(idx_perv,j))
            max_rel_error_liq = max(max_rel_error_liq, diff / max(abs(h2osoi_liq(c,j)), 1.0e-20_r8))
            max_error_vol     = max(max_error_vol,     abs(h2osoi_vol(c,j) - h2osoi_vol_2d(idx_perv,j)))
            if (diff > max_error_liq) then
              max_error_liq    = diff
              max_err_c        = c
              max_err_j        = j
              max_err_idx_perv = idx_perv
              max_err_elm      = h2osoi_liq(c,j)
              max_err_uxx      = h2osoi_liq_2d(idx_perv,j)
            end if
          end do
        end if
      end do

      write(iulog,'(A,I6,A,ES12.4,A,ES12.4,A)') 'Max error in soil water (h2osoi_liq) step=', get_nstep(), &
           ': ', max_error_liq, ' (rel: ', max_rel_error_liq, ')'
      if (max_error_liq > 1.0e-6_r8) then
        write(iulog,*) 'ERROR: Max soil liquid water error exceeds threshold!'
        write(iulog,*) '  ELM column c          :', max_err_c
        write(iulog,*) '  URBANxx idx_perv      :', max_err_idx_perv
        write(iulog,*) '  Layer j               :', max_err_j
        write(iulog,*) '  ELM    h2osoi_liq(c,j):', max_err_elm
        write(iulog,*) '  URBANxx h2osoi_liq    :', max_err_uxx
        write(iulog,*) '  Abs difference        :', max_error_liq
        call exit(0)
      end if

      deallocate(h2osoi_liq_2d)
      deallocate(h2osoi_vol_2d)

    end associate

  end subroutine urbanxx_soilWater_check

end module UrbanxxSoilWaterMod
