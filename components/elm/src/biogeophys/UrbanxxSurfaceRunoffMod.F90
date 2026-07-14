module UrbanxxSurfaceRunoffMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute surface runoff for all five urban surface types via URBANxx.
  ! ELM's SurfaceRunoff logic is replicated inside UrbanComputeSurfaceRunoff
  ! given the assumptions listed in SurfaceRunoff_Porting_Plan.md.
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
  real(c_double) , allocatable, target :: in_qflx_snow_melt_roof(:)
  real(c_double) , allocatable, target :: in_qflx_snow_melt_imperv(:)
  real(c_double) , allocatable, target :: in_qflx_snow_melt_perv(:)
  real(c_double) , allocatable, target :: in_qflx_flood_roof(:)
  real(c_double) , allocatable, target :: in_qflx_flood_imperv(:)
  real(c_double) , allocatable, target :: in_qflx_flood_perv(:)
  real(c_double) , allocatable, target :: in_wtfact_perv(:)
  real(c_double) , allocatable, target :: in_fover_perv(:)
  real(c_double) , allocatable, target :: in_frost_table_perv(:)
  real(c_double) , allocatable, target :: in_zwt_perched_perv(:)

  ! Persistent output buffers (allocated once in init)
  real(c_double) , allocatable, target, public :: out_qflx_surf_roof(:)
  real(c_double) , allocatable, target, public :: out_qflx_surf_imperv(:)
  real(c_double) , allocatable, target, public :: out_qflx_surf_perv(:)
  real(c_double) , allocatable, target, public :: out_qflx_surf_sunwall(:)
  real(c_double) , allocatable, target, public :: out_qflx_surf_shadewall(:)
  real(c_double) , allocatable, target, public :: out_top_h2osoi_liq_roof(:)
  real(c_double) , allocatable, target, public :: out_top_h2osoi_liq_imperv(:)

  public :: urbanxx_surfaceRunoff_init
  public :: urbanxx_surfaceRunoff
  public :: urbanxx_surfaceRunoff_check

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_surfaceRunoff_init(num_urbanl)
    !
    ! !DESCRIPTION:
    ! Allocate persistent buffers for surface runoff computation.
    ! Called once during initialization.
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl

    ! Input buffers
    allocate(in_qflx_snow_melt_roof(num_urbanl))
    allocate(in_qflx_snow_melt_imperv(num_urbanl))
    allocate(in_qflx_snow_melt_perv(num_urbanl))
    allocate(in_qflx_flood_roof(num_urbanl))
    allocate(in_qflx_flood_imperv(num_urbanl))
    allocate(in_qflx_flood_perv(num_urbanl))
    allocate(in_wtfact_perv(num_urbanl))
    allocate(in_fover_perv(num_urbanl))
    allocate(in_frost_table_perv(num_urbanl))
    allocate(in_zwt_perched_perv(num_urbanl))

    ! Output buffers
    allocate(out_qflx_surf_roof(num_urbanl))
    allocate(out_qflx_surf_imperv(num_urbanl))
    allocate(out_qflx_surf_perv(num_urbanl))
    allocate(out_qflx_surf_sunwall(num_urbanl))
    allocate(out_qflx_surf_shadewall(num_urbanl))
    allocate(out_top_h2osoi_liq_roof(num_urbanl))
    allocate(out_top_h2osoi_liq_imperv(num_urbanl))

  end subroutine urbanxx_surfaceRunoff_init

  !-----------------------------------------------------------------------
  subroutine urbanxx_surfaceRunoff(num_urbanl, num_urbanc, filter_urbanc, &
       soilhydrology_vars, soilstate_vars, dtime)
    !
    ! !DESCRIPTION:
    ! Set surface runoff inputs and call UrbanComputeSurfaceRunoff. Outputs
    ! are retrieved into module-level buffers for urbanxx_surfaceRunoff_check
    ! only; never written back into ELM arrays.
    !
    use SoilHydrologyType  , only : soilhydrology_type
    use SoilStateType      , only : soilstate_type
    use ColumnType         , only : col_pp
    use column_varcon      , only : icol_roof, icol_road_imperv, icol_road_perv
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int)          , intent(in) :: num_urbanl
    integer(c_int)          , intent(in) :: num_urbanc
    integer                 , intent(in) :: filter_urbanc(:)
    type(soilhydrology_type), intent(in) :: soilhydrology_vars
    type(soilstate_type)    , intent(in) :: soilstate_vars
    real(r8)                , intent(in) :: dtime
    !
    ! !LOCAL VARIABLES:
    integer(c_int) :: status
    integer        :: fc, c, idx_roof, idx_imperv, idx_perv

    ! Pack ELM's qflx_snomelt (snow melt contribution to qflx_top_soil;
    ! ELM reference: SnowHydrologyMod.F90:472) for roof, impervious road,
    ! and pervious road columns, and set into URBANxx before computing
    ! surface runoff. atm_rain(l) (ForcRain) alone does not capture
    ! melt-water released from a fully-melting initial snowpack.
    idx_roof   = 0
    idx_imperv = 0
    idx_perv   = 0
    do fc = 1, num_urbanc
       c = filter_urbanc(fc)
       if (col_pp%itype(c) == icol_roof) then
          idx_roof = idx_roof + 1
          in_qflx_snow_melt_roof(idx_roof) = col_wf%qflx_snomelt(c)
          in_qflx_flood_roof(idx_roof)     = col_wf%qflx_floodc(c)
       else if (col_pp%itype(c) == icol_road_imperv) then
          idx_imperv = idx_imperv + 1
          in_qflx_snow_melt_imperv(idx_imperv) = col_wf%qflx_snomelt(c)
          in_qflx_flood_imperv(idx_imperv)     = col_wf%qflx_floodc(c)
       else if (col_pp%itype(c) == icol_road_perv) then
          idx_perv = idx_perv + 1
          in_qflx_snow_melt_perv(idx_perv) = col_wf%qflx_snomelt(c)
          in_qflx_flood_perv(idx_perv)     = col_wf%qflx_floodc(c)
          in_wtfact_perv(idx_perv)         = soilstate_vars%wtfact_col(c)
          in_fover_perv(idx_perv)          = soilhydrology_vars%fover(col_pp%gridcell(c))
          in_frost_table_perv(idx_perv)    = soilhydrology_vars%frost_table_col(c)
          in_zwt_perched_perv(idx_perv)    = soilhydrology_vars%zwt_perched_col(c)
       end if
    end do
    call UrbanSetQflxSnowMeltRoof(urbanxx, c_loc(in_qflx_snow_melt_roof), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetQflxSnowMeltImperviousRoad(urbanxx, c_loc(in_qflx_snow_melt_imperv), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetQflxSnowMeltPerviousRoad(urbanxx, c_loc(in_qflx_snow_melt_perv), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Flood water flux (river-routing-model backflow) is added
    ! unconditionally to qflx_surf for every non-wall urban column (ELM
    ! reference: SoilHydrologyMod.F90:247-248).
    call UrbanSetQflxFloodRoof(urbanxx, c_loc(in_qflx_flood_roof), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetQflxFloodImperviousRoad(urbanxx, c_loc(in_qflx_flood_imperv), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetQflxFloodPerviousRoad(urbanxx, c_loc(in_qflx_flood_perv), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Pervious road fsat inputs (ELM reference: SoilHydrologyMod.F90:150-191).
    ! These were previously never set, leaving URBANxx's Wtfact/Fover/
    ! FrostTable/ZwtPerched Kokkos views at their zero-initialized default,
    ! which made fsat (and therefore qflx_surf) identically zero for
    ! pervious road every timestep.
    call UrbanSetWtfactPerviousRoad(urbanxx, c_loc(in_wtfact_perv), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetFoverPerviousRoad(urbanxx, c_loc(in_fover_perv), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetFrostTablePerviousRoad(urbanxx, c_loc(in_frost_table_perv), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanSetZwtPerchedPerviousRoad(urbanxx, c_loc(in_zwt_perched_perv), &
         num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Compute surface runoff
    call UrbanComputeSurfaceRunoff(urbanxx, real(dtime, c_double), status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Retrieve QflxSurf outputs — all five surfaces (diagnostic only; see
    ! urbanxx_surfaceRunoff_check, which compares these against ELM's own
    ! col_wf%qflx_surf/col_ws%h2osoi_liq. Never write URBANxx output back
    ! into ELM data structures.)
    call UrbanGetQflxSurfRoof(urbanxx, c_loc(out_qflx_surf_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxSurfImperviousRoad(urbanxx, c_loc(out_qflx_surf_imperv), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxSurfPerviousRoad(urbanxx, c_loc(out_qflx_surf_perv), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxSurfSunlitWall(urbanxx, c_loc(out_qflx_surf_sunwall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetQflxSurfShadedWall(urbanxx, c_loc(out_qflx_surf_shadewall), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    ! Retrieve updated TopH2OSoiLiq (roof and impervious road)
    call UrbanGetTopH2OSoiLiqRoof(urbanxx, c_loc(out_top_h2osoi_liq_roof), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

    call UrbanGetTopH2OSoiLiqImperviousRoad(urbanxx, c_loc(out_top_h2osoi_liq_imperv), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(iam, __LINE__, status)

  end subroutine urbanxx_surfaceRunoff

  !-----------------------------------------------------------------------
  subroutine urbanxx_surfaceRunoff_check(num_urbanl, num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Compare URBANxx surface runoff outputs against ELM values.
    ! Logs max absolute difference per surface type.
    !
    use ColumnType     , only : col_pp
    use column_varcon  , only : icol_roof, icol_road_imperv, icol_road_perv, &
                                icol_sunwall, icol_shadewall
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer       , intent(in) :: filter_urbanc(:)
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c
    integer  :: idx_roof, idx_imperv, idx_perv, idx_sunwall, idx_shadewall
    real(r8) :: max_err_roof,      max_rel_err_roof
    real(r8) :: max_err_imperv,    max_rel_err_imperv
    real(r8) :: max_err_perv,      max_rel_err_perv
    real(r8) :: max_err_sunwall,   max_rel_err_sunwall
    real(r8) :: max_err_shadewall, max_rel_err_shadewall
    real(r8) :: abs_err, rel_err

    associate( &
         qflx_surf  => col_wf%qflx_surf , &
         h2osoi_liq => col_ws%h2osoi_liq  &
         )

      max_err_roof          = 0._r8 ;  max_rel_err_roof      = 0._r8
      max_err_imperv        = 0._r8 ;  max_rel_err_imperv    = 0._r8
      max_err_perv          = 0._r8 ;  max_rel_err_perv      = 0._r8
      max_err_sunwall       = 0._r8 ;  max_rel_err_sunwall   = 0._r8
      max_err_shadewall     = 0._r8 ;  max_rel_err_shadewall = 0._r8

      idx_roof      = 0
      idx_imperv    = 0
      idx_perv      = 0
      idx_sunwall   = 0
      idx_shadewall = 0
      do fc = 1, num_urbanc
        c = filter_urbanc(fc)
        select case (col_pp%itype(c))
        case (icol_roof)
          idx_roof  = idx_roof + 1
          abs_err   = abs(qflx_surf(c) - out_qflx_surf_roof(idx_roof))
          rel_err   = abs_err / max(abs(qflx_surf(c)), 1.0e-20_r8)
          max_err_roof     = max(max_err_roof,     abs_err)
          max_rel_err_roof = max(max_rel_err_roof, rel_err)
          if (abs_err > 1.0e-10_r8) then
            write(iulog,*) 'ERROR: qflx_surf mismatch (roof) at c=', c
            write(iulog,*) '  ELM=', qflx_surf(c), ' URBANxx=', out_qflx_surf_roof(idx_roof)
            write(iulog,*) '  abs_err=', abs_err, ' rel_err=', rel_err
            call exit(0)
          end if
        case (icol_road_imperv)
          idx_imperv  = idx_imperv + 1
          abs_err     = abs(qflx_surf(c) - out_qflx_surf_imperv(idx_imperv))
          rel_err     = abs_err / max(abs(qflx_surf(c)), 1.0e-20_r8)
          max_err_imperv     = max(max_err_imperv,     abs_err)
          max_rel_err_imperv = max(max_rel_err_imperv, rel_err)
          if (abs_err > 1.0e-10_r8) then
            write(iulog,*) 'ERROR: qflx_surf mismatch (imperv road) at c=', c
            write(iulog,*) '  ELM=', qflx_surf(c), ' URBANxx=', out_qflx_surf_imperv(idx_imperv)
            write(iulog,*) '  abs_err=', abs_err, ' rel_err=', rel_err
            call exit(0)
          end if
        case (icol_road_perv)
          idx_perv  = idx_perv + 1
          abs_err   = abs(qflx_surf(c) - out_qflx_surf_perv(idx_perv))
          rel_err   = abs_err / max(abs(qflx_surf(c)), 1.0e-20_r8)
          max_err_perv     = max(max_err_perv,     abs_err)
          max_rel_err_perv = max(max_rel_err_perv, rel_err)
          if (abs_err > 1.0e-10_r8) then
            write(iulog,*) 'ERROR: qflx_surf mismatch (pervious road) at c=', c
            write(iulog,*) '  ELM=', qflx_surf(c), ' URBANxx=', out_qflx_surf_perv(idx_perv)
            write(iulog,*) '  abs_err=', abs_err, ' rel_err=', rel_err
            call exit(0)
          end if
        case (icol_sunwall)
          idx_sunwall = idx_sunwall + 1
          abs_err     = abs(qflx_surf(c) - out_qflx_surf_sunwall(idx_sunwall))
          rel_err     = abs_err / max(abs(qflx_surf(c)), 1.0e-20_r8)
          max_err_sunwall     = max(max_err_sunwall,     abs_err)
          max_rel_err_sunwall = max(max_rel_err_sunwall, rel_err)
          if (abs_err > 1.0e-10_r8) then
            write(iulog,*) 'ERROR: qflx_surf mismatch (sunlit wall) at c=', c
            write(iulog,*) '  ELM=', qflx_surf(c), ' URBANxx=', out_qflx_surf_sunwall(idx_sunwall)
            write(iulog,*) '  abs_err=', abs_err, ' rel_err=', rel_err
            call exit(0)
          end if
        case (icol_shadewall)
          idx_shadewall = idx_shadewall + 1
          abs_err       = abs(qflx_surf(c) - out_qflx_surf_shadewall(idx_shadewall))
          rel_err       = abs_err / max(abs(qflx_surf(c)), 1.0e-20_r8)
          max_err_shadewall     = max(max_err_shadewall,     abs_err)
          max_rel_err_shadewall = max(max_rel_err_shadewall, rel_err)
          if (abs_err > 1.0e-10_r8) then
            write(iulog,*) 'ERROR: qflx_surf mismatch (shaded wall) at c=', c
            write(iulog,*) '  ELM=', qflx_surf(c), ' URBANxx=', out_qflx_surf_shadewall(idx_shadewall)
            write(iulog,*) '  abs_err=', abs_err, ' rel_err=', rel_err
            call exit(0)
          end if
        end select
      end do

      write(iulog,*) 'Max error in qflx_surf (roof)        : ', max_err_roof,      ' (rel: ', max_rel_err_roof,      ')'
      write(iulog,*) 'Max error in qflx_surf (imperv road) : ', max_err_imperv,    ' (rel: ', max_rel_err_imperv,    ')'
      write(iulog,*) 'Max error in qflx_surf (pervious rd) : ', max_err_perv,      ' (rel: ', max_rel_err_perv,      ')'
      write(iulog,*) 'Max error in qflx_surf (sunlit wall) : ', max_err_sunwall,   ' (rel: ', max_rel_err_sunwall,   ')'
      write(iulog,*) 'Max error in qflx_surf (shaded wall) : ', max_err_shadewall, ' (rel: ', max_rel_err_shadewall, ')'

    end associate

  end subroutine urbanxx_surfaceRunoff_check

end module UrbanxxSurfaceRunoffMod
