module UrbanxxDiagnosticControlMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Namelist-driven control for URBANxx's opt-in validation modes. Several
  ! ELM-vs-URBANxx physics-parity checks compare independently-computed
  ! values to a tight (1e-6 to 1e-9, depending on quantity) tolerance and
  ! abort the run on mismatch. One specific, thoroughly-investigated
  ! residual (pervious-road column near a diurnal freeze/thaw boundary,
  ! first appears ~nstep=203 in the modified_f19_f19 case) is understood to
  ! be a genuine, unavoidable consequence of Fortran and C++ not
  ! reproducing bit-identical rounding for algebraically-identical
  ! floating-point formulas — literally a single-ULP difference that
  ! compounds over repeated freeze/thaw events — not a coding defect. See
  ! plans/ELM_URBANxx_Known_Gaps_And_Future_Work.md, "Round 1" through
  ! "Round 6" for the full investigation and root-cause trace.
  !
  ! Two independent, namelist-controlled modes let a run continue past
  ! this kind of tolerance exceedance instead of aborting:
  !
  !   urbanxx_softfail (logical, default .false.): bounded soft-fail.
  !     Logs a warning and continues past an exceedance, but still aborts
  !     if the exceedance is large enough to plausibly be a *different,
  !     genuine* problem, governed by a ceiling expressed as a multiplier
  !     on each check's own tolerance (urbanxx_softfail_ceiling_factor,
  !     default 1.0e6) so it scales naturally per quantity/unit without
  !     needing a separate absolute ceiling per check.
  !
  !   urbanxx_always_continue (logical, default .false.): unconditional
  !     continue. Every check logs its error and the run always continues
  !     to completion, regardless of exceedance magnitude. Intended for
  !     full-length validation sweeps that want to characterize every
  !     divergence in one run rather than stopping at the first one.
  !
  ! Default is fully disabled (strict/abort, matching every check's
  ! original behavior) — this module has zero effect on production runs
  ! unless explicitly opted into via user_nl_elm.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use elm_varctl,   only : urbanxx_softfail, urbanxx_softfail_ceiling_factor, &
                           urbanxx_always_continue

  implicit none

  private

  public :: urbanxx_softfail_permits

contains

  !-----------------------------------------------------------------------
  logical function urbanxx_softfail_permits(error, tol) result(permits)
    !
    ! !DESCRIPTION:
    ! Convenience wrapper for check sites: returns .true. if the run
    ! should log a warning and continue rather than abort for this
    ! exceedance — either because urbanxx_always_continue is set (always
    ! permits), or because urbanxx_softfail is set and this exceedance is
    ! within the tolerance ceiling. Note IEEE NaN comparisons are always
    ! .false., so a NaN error already fails the caller's own
    ! "error > tol" trip condition before this function is even reached —
    ! this is a pre-existing property of every check site, not something
    ! either mode changes.
    !
    real(r8), intent(in) :: error
    real(r8), intent(in) :: tol

    permits = urbanxx_always_continue .or. &
         (urbanxx_softfail .and. (error <= tol * urbanxx_softfail_ceiling_factor))

  end function urbanxx_softfail_permits

end module UrbanxxDiagnosticControlMod
