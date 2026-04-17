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
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use spmdMod              , only : masterproc, iam
  use elm_varctl           , only : iulog
  use ColumnDataType       , only : col_wf
  use UrbanxxInstanceMod   , only : urbanxx

  implicit none

  private

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
    implicit none
    integer(c_int), intent(in) :: num_urbanl

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
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer       , intent(in) :: filter_urbanc(:)           ! urban column filter
    !
    ! !LOCAL VARIABLES:
    integer(c_int)                       :: status

    call UrbanComputeInfiltration(urbanxx, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetInfiltrationFluxPerviousRoad(urbanxx, c_loc(out_qflxInfl), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

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
