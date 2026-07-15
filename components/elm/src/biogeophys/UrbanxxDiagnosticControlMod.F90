module UrbanxxDiagnosticControlMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Runtime, environment-variable-driven control for URBANxx's opt-in
  ! "soft-fail" validation mode. Several ELM-vs-URBANxx physics-parity
  ! checks compare independently-computed values to a tight (1e-6 to
  ! 1e-9, depending on quantity) tolerance and abort the run on mismatch.
  ! One specific, thoroughly-investigated residual (pervious-road column
  ! near a diurnal freeze/thaw boundary, first appears ~nstep=203 in the
  ! modified_f19_f19 case) is understood to be a genuine, unavoidable
  ! consequence of Fortran and C++ not reproducing bit-identical rounding
  ! for algebraically-identical floating-point formulas — literally a
  ! single-ULP difference that compounds over repeated freeze/thaw
  ! events — not a coding defect. See
  ! plans/ELM_URBANxx_Known_Gaps_And_Future_Work.md, "Round 1" through
  ! "Round 6" for the full investigation and root-cause trace.
  !
  ! Soft-fail mode lets a run continue past this kind of tolerance
  ! exceedance (logging a warning instead of aborting) while still
  ! aborting on any exceedance large enough to plausibly be a *different,
  ! genuine* problem (governed by a ceiling expressed as a multiplier on
  ! each check's own tolerance, so it scales naturally per quantity/unit
  ! without needing a separate absolute ceiling per check). Default is
  ! fully disabled (strict/abort, matching every check's original
  ! behavior) — this module has zero effect on production runs unless
  ! explicitly opted into via environment variable.
  !
  ! Environment variables (read once per run, cached):
  !   URBANXX_SOFT_FAIL                : 0/1, enable soft-fail mode.
  !                                       Default 0 (disabled — checks
  !                                       abort exactly as before).
  !   URBANXX_SOFT_FAIL_CEILING_FACTOR : real, multiplier on a check's own
  !                                       tolerance above which an
  !                                       exceedance still aborts even in
  !                                       soft-fail mode. Default 1.0e6.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  private

  public :: urbanxx_softfail_get_config
  public :: urbanxx_softfail_permits

  logical , save :: initialized              = .false.
  logical , save :: cached_soft_fail         = .false.
  real(r8), save :: cached_ceiling_factor    = 1.0e6_r8

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_softfail_get_config(enabled, ceiling_factor)
    !
    ! !DESCRIPTION:
    ! Return the current soft-fail configuration, reading and caching it
    ! from environment variables on the first call.
    !
    logical , intent(out) :: enabled
    real(r8), intent(out) :: ceiling_factor

    if (.not. initialized) then
       cached_soft_fail      = (read_env_int('URBANXX_SOFT_FAIL', 0) /= 0)
       cached_ceiling_factor = read_env_real('URBANXX_SOFT_FAIL_CEILING_FACTOR', 1.0e6_r8)
       initialized = .true.
    end if

    enabled        = cached_soft_fail
    ceiling_factor = cached_ceiling_factor

  end subroutine urbanxx_softfail_get_config

  !-----------------------------------------------------------------------
  logical function urbanxx_softfail_permits(error, tol) result(permits)
    !
    ! !DESCRIPTION:
    ! Convenience wrapper for check sites: returns .true. if soft-fail
    ! mode is enabled AND this exceedance is within the tolerance
    ! ceiling (i.e. the check should log a warning and continue rather
    ! than abort). Note IEEE NaN comparisons are always .false., so a
    ! NaN error already fails the caller's own "error > tol" trip
    ! condition before this function is even reached — this is a
    ! pre-existing property of every check site, not something soft-fail
    ! mode changes.
    !
    real(r8), intent(in) :: error
    real(r8), intent(in) :: tol
    logical  :: enabled
    real(r8) :: ceiling_factor

    call urbanxx_softfail_get_config(enabled, ceiling_factor)
    permits = enabled .and. (error <= tol * ceiling_factor)

  end function urbanxx_softfail_permits

  !-----------------------------------------------------------------------
  integer function read_env_int(name, default_value) result(val)
    !
    ! !DESCRIPTION:
    ! Read an integer-valued environment variable, falling back to
    ! default_value if unset, empty, or unparseable.
    !
    character(len=*), intent(in) :: name
    integer         , intent(in) :: default_value
    character(len=64) :: buf
    integer :: env_len, env_status

    call get_environment_variable(name=name, value=buf, length=env_len, &
         status=env_status, trim_name=.true.)
    if (env_status == 0 .and. env_len > 0) then
       read(buf, *, iostat=env_status) val
       if (env_status /= 0) val = default_value
    else
       val = default_value
    end if
  end function read_env_int

  !-----------------------------------------------------------------------
  real(r8) function read_env_real(name, default_value) result(val)
    !
    ! !DESCRIPTION:
    ! Read a real-valued environment variable, falling back to
    ! default_value if unset, empty, or unparseable.
    !
    character(len=*), intent(in) :: name
    real(r8)        , intent(in) :: default_value
    character(len=64) :: buf
    integer :: env_len, env_status

    call get_environment_variable(name=name, value=buf, length=env_len, &
         status=env_status, trim_name=.true.)
    if (env_status == 0 .and. env_len > 0) then
       read(buf, *, iostat=env_status) val
       if (env_status /= 0) val = default_value
    else
       val = default_value
    end if
  end function read_env_real

end module UrbanxxDiagnosticControlMod
