module UrbanxxInfiltrationMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute infiltration flux for pervious road urban columns in the
  ! URBANxx model.  The simplified urban formula (snl >= 0 always true for
  ! urban) is:
  !   qflx_infl = qflx_top_soil - qflx_surf - qflx_evap_grnd
  ! This module sets the required inputs, calls UrbanComputeInfiltration,
  ! and retrieves the resulting QflxInfl into the URBANxx pervious-road
  ! data structure so that UrbanComputeHydrology can consume it.
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
  real(c_double) , allocatable, target :: h2oLiq(:)        ! h2osoi_liq, flattened (num_urbanl x nlevgrnd)
  real(c_double) , allocatable, target :: h2oIce(:)        ! h2osoi_ice, flattened (num_urbanl x nlevgrnd)

  ! Persistent output buffer (allocated once in init)
  real(c_double) , allocatable, target, public :: out_qflxInfl(:)  ! output qflx_infl from URBANxx (num_urbanl)

  public :: urbanxx_infiltration_init
  public :: urbanxx_infiltration
  public :: urbanxx_infiltration_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_infiltration_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for infiltration computation.
    ! Called once during initialization.
    !
    use elm_varpar, only : nlevgrnd
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int) :: totalSize

    totalSize = num_urbanl * nlevgrnd

    ! Input buffers — 2D flattened
    allocate(h2oLiq(totalSize))
    allocate(h2oIce(totalSize))

    ! Output buffer — 1D
    allocate(out_qflxInfl(num_urbanl))

  end subroutine urbanxx_infiltration_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_infiltration(num_urbanl, num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Set infiltration inputs for pervious road and compute infiltration
    ! flux via URBANxx.
    !
    use ColumnType           , only : col_pp
    use column_varcon        , only : icol_road_perv
    use elm_varpar           , only : nlevgrnd
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer       , intent(in) :: filter_urbanc(:)           ! urban column filter
    !
    ! !LOCAL VARIABLES:
    integer(c_int)                       :: status
    integer                              :: fc, c, j, idx, idx_perv, nlevbed
    integer(c_int)                       :: totalSize
    integer(c_int), dimension(2)         :: size2D
    logical(c_bool)                      :: isLayoutLeft

    associate(                               &
         qflx_surf      => col_wf%qflx_surf      , & ! Input: [real(r8) (:)] surface runoff (mm H2O /s)
         qflx_evap_grnd => col_wf%qflx_evap_grnd , & ! Input: [real(r8) (:)] ground surface evaporation rate (mm H2O/s)
         h2osoi_liq     => col_ws%h2osoi_liq     , & ! Input: [real(r8) (:,:)] liquid water (kg/m2)
         h2osoi_ice     => col_ws%h2osoi_ice     , & ! Input: [real(r8) (:,:)] ice lens (kg/m2)
         nlev2bed       => col_pp%nlevbed           & ! Input: [integer (:)] number of layers to bedrock
         )

      ! Pack 2D soil water buffers for pervious road columns
      totalSize = num_urbanl * nlevgrnd
      size2D(1) = num_urbanl
      size2D(2) = nlevgrnd

      isLayoutLeft = UrbanKokkosIsLayoutLeft()

      if (isLayoutLeft) then
        ! LayoutLeft: first dimension (landunits) varies fastest
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
              else
                h2oLiq(idx) = 0.0_r8
                h2oIce(idx) = 0.0_r8
              end if
            end if
          end do
        end do
      else
        ! LayoutRight: last dimension (layers) varies fastest
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
              else
                h2oLiq(idx) = 0.0_r8
                h2oIce(idx) = 0.0_r8
              end if
            end do
          end if
        end do
      end if

      call UrbanSetSoilLiquidWaterForPerviousRoad(urbanxx, c_loc(h2oLiq), &
           size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanSetSoilIceContentForPerviousRoad(urbanxx, c_loc(h2oIce), &
           size2D, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanComputeInfiltration(urbanxx, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

      call UrbanGetInfiltrationFluxPerviousRoad(urbanxx, c_loc(out_qflxInfl), &
           num_urbanl, status)
      if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    end associate

  end subroutine urbanxx_infiltration

  !-----------------------------------------------------------------------
  subroutine urbanxx_infiltration_check(num_urbanl, num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Compare URBANxx infiltration output against ELM values for pervious
    ! road columns.  Reports max absolute and relative differences.
    !
    use ColumnType     , only : col_pp
    use column_varcon  , only : icol_road_perv
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer       , intent(in) :: filter_urbanc(:)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, idx_perv
    real(r8) :: max_abs_err, max_rel_err, abs_err, rel_err

    associate( &
         qflx_infl => col_wf%qflx_infl &
         )

      max_abs_err = 0._r8
      max_rel_err = 0._r8
      idx_perv    = 0

      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          abs_err = abs(qflx_infl(c) - out_qflxInfl(idx_perv))
          rel_err = abs_err / max(abs(qflx_infl(c)), 1.0e-20_r8)
          max_abs_err = max(max_abs_err, abs_err)
          max_rel_err = max(max_rel_err, rel_err)
          if (max_abs_err > 1.0e-9_r8) then
            write(iulog,*)'Max error in qflx_infl: ', max_abs_err, ' at column ', c, &
                 ' itype = ', col_pp%itype(c), ' exceeds tolerance'
            call exit(0)
          end if
        end if
      end do

      write(iulog,*) 'Max error in infiltration flux      : ', max_abs_err, ' (rel: ', max_rel_err, ')'

    end associate

  end subroutine urbanxx_infiltration_check

end module UrbanxxInfiltrationMod
