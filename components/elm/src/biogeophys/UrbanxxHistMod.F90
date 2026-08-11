module UrbanxxHistMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! History output for URBANxx state and fluxes.  All fields carry the
  ! _XX suffix to distinguish them from ELM's own solver output.
  !
  ! Controlled by the namelist flag urbanxx_output (default .false.).
  ! When .false. no arrays are allocated, no hist_addfld calls are made,
  ! and urbanxx_hist_update is a no-op.
  !
  ! urbanxx_hist_init must be called from initialize2 (after
  ! elm_inst_biogeophys but before hist_htapes_build) so that hist_addfld
  ! registrations happen before the tape field list is finalized.
  !-----------------------------------------------------------------------

  use iso_c_binding
  use urban_mod
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use elm_varcon    , only : spval
  use elm_varctl    , only : urbanxx_output, use_urbanxx
  use elm_varpar    , only : nlevurb, nlevgrnd, nlevsno
  use histFileMod   , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use decompMod     , only : bounds_type
  use UrbanxxInstanceMod, only : urbanxx
  use urban_kokkos_interface, only : UrbanKokkosIsLayoutLeft

  implicit none
  private

  public :: urbanxx_hist_init
  public :: urbanxx_hist_update
  public :: urbanxx_hist_update_radiation

  ! ---- module-level output arrays (allocated only when urbanxx_output=.true.) ----

  ! 2D layer temperature: (begc:endc, 1:nlevurb)
  real(r8), allocatable, target :: out_tsoi_roof(:,:)
  real(r8), allocatable, target :: out_tsoi_sunwall(:,:)
  real(r8), allocatable, target :: out_tsoi_shadwall(:,:)
  ! 2D layer temperature: (begc:endc, 1:nlevgrnd)
  real(r8), allocatable, target :: out_tsoi_improad(:,:)
  real(r8), allocatable, target :: out_tsoi_perroad(:,:)

  ! 2D snow layer fields: (begc:endc, 1:nlevsno)
  real(r8), allocatable, target :: out_sno_t_roof(:,:)
  real(r8), allocatable, target :: out_sno_t_improad(:,:)
  real(r8), allocatable, target :: out_sno_t_perroad(:,:)
  real(r8), allocatable, target :: out_sno_liq_roof(:,:)
  real(r8), allocatable, target :: out_sno_liq_improad(:,:)
  real(r8), allocatable, target :: out_sno_liq_perroad(:,:)
  real(r8), allocatable, target :: out_sno_ice_roof(:,:)
  real(r8), allocatable, target :: out_sno_ice_improad(:,:)
  real(r8), allocatable, target :: out_sno_ice_perroad(:,:)
  real(r8), allocatable, target :: out_sno_dz_roof(:,:)
  real(r8), allocatable, target :: out_sno_dz_improad(:,:)
  real(r8), allocatable, target :: out_sno_dz_perroad(:,:)

  ! 2D pervious road soil water: (begc:endc, 1:nlevgrnd)
  real(r8), allocatable, target :: out_soilliq_perroad(:,:)
  real(r8), allocatable, target :: out_soilice_perroad(:,:)

  ! 1D snow state scalars: (begc:endc)
  real(r8), allocatable, target :: out_h2osno_roof(:)
  real(r8), allocatable, target :: out_h2osno_improad(:)
  real(r8), allocatable, target :: out_snol_roof(:)
  real(r8), allocatable, target :: out_snol_improad(:)
  real(r8), allocatable, target :: out_snol_perroad(:)

  ! 1D building/canyon scalars assigned to roof column: (begc:endc)
  real(r8), allocatable, target :: out_t_canyon(:)
  real(r8), allocatable, target :: out_t_bld(:)

  ! Phase 2 — Flux and energy balance fields: (begc:endc)
  real(r8), allocatable, target :: out_eflx_sg_roof(:)
  real(r8), allocatable, target :: out_eflx_sg_sunwall(:)
  real(r8), allocatable, target :: out_eflx_sg_shadwall(:)
  real(r8), allocatable, target :: out_eflx_sg_improad(:)
  real(r8), allocatable, target :: out_eflx_sg_perroad(:)
  real(r8), allocatable, target :: out_xmf_roof(:)
  real(r8), allocatable, target :: out_xmf_improad(:)
  real(r8), allocatable, target :: out_xmf_perroad(:)
  real(r8), allocatable, target :: out_netlw_roof(:)
  real(r8), allocatable, target :: out_netlw_sunwall(:)
  real(r8), allocatable, target :: out_netlw_shadwall(:)
  real(r8), allocatable, target :: out_netlw_improad(:)
  real(r8), allocatable, target :: out_netlw_perroad(:)
  real(r8), allocatable, target :: out_fsh_roof(:)
  real(r8), allocatable, target :: out_fsh_sunwall(:)
  real(r8), allocatable, target :: out_fsh_shadwall(:)
  real(r8), allocatable, target :: out_fsh_improad(:)
  real(r8), allocatable, target :: out_fsh_perroad(:)
  real(r8), allocatable, target :: out_qflxsnomelt_roof(:)
  real(r8), allocatable, target :: out_qflxsnomelt_improad(:)
  real(r8), allocatable, target :: out_qflxsnomelt_perroad(:)
  real(r8), allocatable, target :: out_qflxsubsnow_roof(:)
  real(r8), allocatable, target :: out_qflxsubsnow_improad(:)
  real(r8), allocatable, target :: out_qflxsubsnow_perroad(:)
  real(r8), allocatable, target :: out_qflxsurf_roof(:)
  real(r8), allocatable, target :: out_qflxsurf_improad(:)
  real(r8), allocatable, target :: out_qflxsurf_perroad(:)

  ! Phase 3 — Net shortwave radiation fields: (begc:endc)
  real(r8), allocatable, target :: out_netsw_roof(:)
  real(r8), allocatable, target :: out_netsw_sunwall(:)
  real(r8), allocatable, target :: out_netsw_shadwall(:)
  real(r8), allocatable, target :: out_netsw_improad(:)
  real(r8), allocatable, target :: out_netsw_perroad(:)

contains

  !-----------------------------------------------------------------------
  subroutine urbanxx_hist_init(bounds)
    !
    ! !DESCRIPTION:
    ! Allocate output arrays and register all URBANxx history fields.
    ! Must be called from initialize2, after elm_inst_biogeophys but before
    ! hist_htapes_build.
    !
    use abortutils, only : endrun
    !
    implicit none
    type(bounds_type), intent(in) :: bounds

    ! Local pointer variables for hist_addfld registration
    real(r8), pointer :: ptr1d(:)
    real(r8), pointer :: ptr2d(:,:)

    integer :: begc, endc

    if (.not. urbanxx_output) return

    if (.not. use_urbanxx) then
       call endrun(msg='urbanxx_hist_init: urbanxx_output=.true. requires use_urbanxx=.true.')
    end if

    begc = bounds%begc
    endc = bounds%endc

    ! --- Allocate and initialize all output arrays ---

    allocate(out_tsoi_roof    (begc:endc, 1:nlevurb));  out_tsoi_roof     = spval
    allocate(out_tsoi_sunwall (begc:endc, 1:nlevurb));  out_tsoi_sunwall  = spval
    allocate(out_tsoi_shadwall(begc:endc, 1:nlevurb));  out_tsoi_shadwall = spval
    allocate(out_tsoi_improad (begc:endc, 1:nlevgrnd)); out_tsoi_improad  = spval
    allocate(out_tsoi_perroad (begc:endc, 1:nlevgrnd)); out_tsoi_perroad  = spval

    allocate(out_sno_t_roof    (begc:endc, 1:nlevsno)); out_sno_t_roof     = spval
    allocate(out_sno_t_improad (begc:endc, 1:nlevsno)); out_sno_t_improad  = spval
    allocate(out_sno_t_perroad (begc:endc, 1:nlevsno)); out_sno_t_perroad  = spval
    allocate(out_sno_liq_roof  (begc:endc, 1:nlevsno)); out_sno_liq_roof   = spval
    allocate(out_sno_liq_improad(begc:endc, 1:nlevsno));out_sno_liq_improad= spval
    allocate(out_sno_liq_perroad(begc:endc, 1:nlevsno));out_sno_liq_perroad= spval
    allocate(out_sno_ice_roof  (begc:endc, 1:nlevsno)); out_sno_ice_roof   = spval
    allocate(out_sno_ice_improad(begc:endc, 1:nlevsno));out_sno_ice_improad= spval
    allocate(out_sno_ice_perroad(begc:endc, 1:nlevsno));out_sno_ice_perroad= spval
    allocate(out_sno_dz_roof   (begc:endc, 1:nlevsno)); out_sno_dz_roof    = spval
    allocate(out_sno_dz_improad(begc:endc, 1:nlevsno)); out_sno_dz_improad = spval
    allocate(out_sno_dz_perroad(begc:endc, 1:nlevsno)); out_sno_dz_perroad = spval

    allocate(out_soilliq_perroad(begc:endc, 1:nlevgrnd)); out_soilliq_perroad = spval
    allocate(out_soilice_perroad(begc:endc, 1:nlevgrnd)); out_soilice_perroad = spval

    allocate(out_h2osno_roof  (begc:endc)); out_h2osno_roof   = spval
    allocate(out_h2osno_improad(begc:endc));out_h2osno_improad = spval
    allocate(out_snol_roof    (begc:endc)); out_snol_roof     = spval
    allocate(out_snol_improad (begc:endc)); out_snol_improad  = spval
    allocate(out_snol_perroad (begc:endc)); out_snol_perroad  = spval

    allocate(out_t_canyon(begc:endc)); out_t_canyon = spval
    allocate(out_t_bld   (begc:endc)); out_t_bld    = spval

    ! Phase 2 allocations
    allocate(out_eflx_sg_roof    (begc:endc)); out_eflx_sg_roof     = spval
    allocate(out_eflx_sg_sunwall (begc:endc)); out_eflx_sg_sunwall  = spval
    allocate(out_eflx_sg_shadwall(begc:endc)); out_eflx_sg_shadwall = spval
    allocate(out_eflx_sg_improad (begc:endc)); out_eflx_sg_improad  = spval
    allocate(out_eflx_sg_perroad (begc:endc)); out_eflx_sg_perroad  = spval
    allocate(out_xmf_roof        (begc:endc)); out_xmf_roof         = spval
    allocate(out_xmf_improad     (begc:endc)); out_xmf_improad      = spval
    allocate(out_xmf_perroad     (begc:endc)); out_xmf_perroad      = spval
    allocate(out_netlw_roof      (begc:endc)); out_netlw_roof        = spval
    allocate(out_netlw_sunwall   (begc:endc)); out_netlw_sunwall     = spval
    allocate(out_netlw_shadwall  (begc:endc)); out_netlw_shadwall    = spval
    allocate(out_netlw_improad   (begc:endc)); out_netlw_improad     = spval
    allocate(out_netlw_perroad   (begc:endc)); out_netlw_perroad     = spval
    allocate(out_fsh_roof        (begc:endc)); out_fsh_roof          = spval
    allocate(out_fsh_sunwall     (begc:endc)); out_fsh_sunwall       = spval
    allocate(out_fsh_shadwall    (begc:endc)); out_fsh_shadwall      = spval
    allocate(out_fsh_improad     (begc:endc)); out_fsh_improad       = spval
    allocate(out_fsh_perroad     (begc:endc)); out_fsh_perroad       = spval
    allocate(out_qflxsnomelt_roof   (begc:endc)); out_qflxsnomelt_roof    = spval
    allocate(out_qflxsnomelt_improad(begc:endc)); out_qflxsnomelt_improad = spval
    allocate(out_qflxsnomelt_perroad(begc:endc)); out_qflxsnomelt_perroad = spval
    allocate(out_qflxsubsnow_roof   (begc:endc)); out_qflxsubsnow_roof    = spval
    allocate(out_qflxsubsnow_improad(begc:endc)); out_qflxsubsnow_improad = spval
    allocate(out_qflxsubsnow_perroad(begc:endc)); out_qflxsubsnow_perroad = spval
    allocate(out_qflxsurf_roof      (begc:endc)); out_qflxsurf_roof       = spval
    allocate(out_qflxsurf_improad   (begc:endc)); out_qflxsurf_improad    = spval
    allocate(out_qflxsurf_perroad   (begc:endc)); out_qflxsurf_perroad    = spval
    ! Phase 3 allocations
    allocate(out_netsw_roof    (begc:endc)); out_netsw_roof     = spval
    allocate(out_netsw_sunwall (begc:endc)); out_netsw_sunwall  = spval
    allocate(out_netsw_shadwall(begc:endc)); out_netsw_shadwall = spval
    allocate(out_netsw_improad (begc:endc)); out_netsw_improad  = spval
    allocate(out_netsw_perroad (begc:endc)); out_netsw_perroad  = spval

    ! --- Register 2D layer temperature fields ---

    ptr2d => out_tsoi_roof
    call hist_addfld2d(fname='TSOI_ROOF_XX', type2d='levurb', units='K', &
         avgflag='A', long_name='URBANxx roof layer temperatures', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr2d => out_tsoi_sunwall
    call hist_addfld2d(fname='TSOI_SUNWALL_XX', type2d='levurb', units='K', &
         avgflag='A', long_name='URBANxx sunlit wall layer temperatures', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr2d => out_tsoi_shadwall
    call hist_addfld2d(fname='TSOI_SHADWALL_XX', type2d='levurb', units='K', &
         avgflag='A', long_name='URBANxx shaded wall layer temperatures', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr2d => out_tsoi_improad
    call hist_addfld2d(fname='TSOI_IMPROAD_XX', type2d='levgrnd', units='K', &
         avgflag='A', long_name='URBANxx impervious road layer temperatures', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr2d => out_tsoi_perroad
    call hist_addfld2d(fname='TSOI_PERROAD_XX', type2d='levgrnd', units='K', &
         avgflag='A', long_name='URBANxx pervious road layer temperatures', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ! --- Register 1D snow state scalars ---

    ptr1d => out_h2osno_roof
    call hist_addfld1d(fname='H2OSNO_ROOF_XX', units='kg/m2', &
         avgflag='A', long_name='URBANxx roof snow water equivalent', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_h2osno_improad
    call hist_addfld1d(fname='H2OSNO_IMPROAD_XX', units='kg/m2', &
         avgflag='A', long_name='URBANxx impervious road snow water equivalent', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_snol_roof
    call hist_addfld1d(fname='SNOL_ROOF_XX', units='1', &
         avgflag='I', long_name='URBANxx roof number of active snow layers', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_snol_improad
    call hist_addfld1d(fname='SNOL_IMPROAD_XX', units='1', &
         avgflag='I', long_name='URBANxx impervious road number of active snow layers', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_snol_perroad
    call hist_addfld1d(fname='SNOL_PERROAD_XX', units='1', &
         avgflag='I', long_name='URBANxx pervious road number of active snow layers', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ! --- Register 2D snow layer fields ---

    ptr2d => out_sno_t_roof
    call hist_addfld2d(fname='SNO_T_ROOF_XX', type2d='levsno', units='K', &
         avgflag='A', long_name='URBANxx roof snow layer temperatures', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_t_improad
    call hist_addfld2d(fname='SNO_T_IMPROAD_XX', type2d='levsno', units='K', &
         avgflag='A', long_name='URBANxx impervious road snow layer temperatures', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_t_perroad
    call hist_addfld2d(fname='SNO_T_PERROAD_XX', type2d='levsno', units='K', &
         avgflag='A', long_name='URBANxx pervious road snow layer temperatures', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_liq_roof
    call hist_addfld2d(fname='SNO_LIQ_ROOF_XX', type2d='levsno', units='kg/m2', &
         avgflag='A', long_name='URBANxx roof snow layer liquid water', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_liq_improad
    call hist_addfld2d(fname='SNO_LIQ_IMPROAD_XX', type2d='levsno', units='kg/m2', &
         avgflag='A', long_name='URBANxx impervious road snow layer liquid water', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_liq_perroad
    call hist_addfld2d(fname='SNO_LIQ_PERROAD_XX', type2d='levsno', units='kg/m2', &
         avgflag='A', long_name='URBANxx pervious road snow layer liquid water', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_ice_roof
    call hist_addfld2d(fname='SNO_ICE_ROOF_XX', type2d='levsno', units='kg/m2', &
         avgflag='A', long_name='URBANxx roof snow layer ice content', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_ice_improad
    call hist_addfld2d(fname='SNO_ICE_IMPROAD_XX', type2d='levsno', units='kg/m2', &
         avgflag='A', long_name='URBANxx impervious road snow layer ice content', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_ice_perroad
    call hist_addfld2d(fname='SNO_ICE_PERROAD_XX', type2d='levsno', units='kg/m2', &
         avgflag='A', long_name='URBANxx pervious road snow layer ice content', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_dz_roof
    call hist_addfld2d(fname='SNO_DZ_ROOF_XX', type2d='levsno', units='m', &
         avgflag='A', long_name='URBANxx roof snow layer thickness', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_dz_improad
    call hist_addfld2d(fname='SNO_DZ_IMPROAD_XX', type2d='levsno', units='m', &
         avgflag='A', long_name='URBANxx impervious road snow layer thickness', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ptr2d => out_sno_dz_perroad
    call hist_addfld2d(fname='SNO_DZ_PERROAD_XX', type2d='levsno', units='m', &
         avgflag='A', long_name='URBANxx pervious road snow layer thickness', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', &
         no_snow_behavior=no_snow_normal, default='inactive')

    ! --- Register 2D pervious road soil water fields ---

    ptr2d => out_soilliq_perroad
    call hist_addfld2d(fname='SOILLIQ_PERROAD_XX', type2d='levgrnd', units='kg/m2', &
         avgflag='A', long_name='URBANxx pervious road soil liquid water', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', default='inactive')

    ptr2d => out_soilice_perroad
    call hist_addfld2d(fname='SOILICE_PERROAD_XX', type2d='levgrnd', units='kg/m2', &
         avgflag='A', long_name='URBANxx pervious road soil ice content', &
         ptr_col=ptr2d, set_nourb=0._r8, c2l_scale_type='urbanf', default='inactive')

    ! --- Register 1D building/canyon scalars (assigned to roof column) ---

    ptr1d => out_t_canyon
    call hist_addfld1d(fname='T_CANYON_XX', units='K', &
         avgflag='A', long_name='URBANxx canyon air temperature', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_t_bld
    call hist_addfld1d(fname='T_BLD_XX', units='K', &
         avgflag='A', long_name='URBANxx building interior temperature', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ! --- Phase 2: Register flux and energy balance fields ---

    ptr1d => out_eflx_sg_roof
    call hist_addfld1d(fname='EFLX_SG_ROOF_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx roof ground heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_eflx_sg_sunwall
    call hist_addfld1d(fname='EFLX_SG_SUNWALL_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx sunlit wall ground heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_eflx_sg_shadwall
    call hist_addfld1d(fname='EFLX_SG_SHADWALL_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx shaded wall ground heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_eflx_sg_improad
    call hist_addfld1d(fname='EFLX_SG_IMPROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx impervious road ground heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_eflx_sg_perroad
    call hist_addfld1d(fname='EFLX_SG_PERROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx pervious road ground heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_xmf_roof
    call hist_addfld1d(fname='XMF_ROOF_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx roof phase change heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_xmf_improad
    call hist_addfld1d(fname='XMF_IMPROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx impervious road phase change heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_xmf_perroad
    call hist_addfld1d(fname='XMF_PERROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx pervious road phase change heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_netlw_roof
    call hist_addfld1d(fname='NETLW_ROOF_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx roof net longwave radiation', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_netlw_sunwall
    call hist_addfld1d(fname='NETLW_SUNWALL_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx sunlit wall net longwave radiation', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_netlw_shadwall
    call hist_addfld1d(fname='NETLW_SHADWALL_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx shaded wall net longwave radiation', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_netlw_improad
    call hist_addfld1d(fname='NETLW_IMPROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx impervious road net longwave radiation', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_netlw_perroad
    call hist_addfld1d(fname='NETLW_PERROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx pervious road net longwave radiation', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_fsh_roof
    call hist_addfld1d(fname='FSH_ROOF_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx roof sensible heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_fsh_sunwall
    call hist_addfld1d(fname='FSH_SUNWALL_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx sunlit wall sensible heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_fsh_shadwall
    call hist_addfld1d(fname='FSH_SHADWALL_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx shaded wall sensible heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_fsh_improad
    call hist_addfld1d(fname='FSH_IMPROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx impervious road sensible heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_fsh_perroad
    call hist_addfld1d(fname='FSH_PERROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx pervious road sensible heat flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_qflxsnomelt_roof
    call hist_addfld1d(fname='QFLX_SNOMELT_ROOF_XX', units='mm/s', &
         avgflag='A', long_name='URBANxx roof snow melt flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_qflxsnomelt_improad
    call hist_addfld1d(fname='QFLX_SNOMELT_IMPROAD_XX', units='mm/s', &
         avgflag='A', long_name='URBANxx impervious road snow melt flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_qflxsnomelt_perroad
    call hist_addfld1d(fname='QFLX_SNOMELT_PERROAD_XX', units='mm/s', &
         avgflag='A', long_name='URBANxx pervious road snow melt flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_qflxsubsnow_roof
    call hist_addfld1d(fname='QFLX_SUBSNOW_ROOF_XX', units='mm/s', &
         avgflag='A', long_name='URBANxx roof snow sublimation flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_qflxsubsnow_improad
    call hist_addfld1d(fname='QFLX_SUBSNOW_IMPROAD_XX', units='mm/s', &
         avgflag='A', long_name='URBANxx impervious road snow sublimation flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_qflxsubsnow_perroad
    call hist_addfld1d(fname='QFLX_SUBSNOW_PERROAD_XX', units='mm/s', &
         avgflag='A', long_name='URBANxx pervious road snow sublimation flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_qflxsurf_roof
    call hist_addfld1d(fname='QFLX_SURF_ROOF_XX', units='mm/s', &
         avgflag='A', long_name='URBANxx roof surface runoff flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_qflxsurf_improad
    call hist_addfld1d(fname='QFLX_SURF_IMPROAD_XX', units='mm/s', &
         avgflag='A', long_name='URBANxx impervious road surface runoff flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_qflxsurf_perroad
    call hist_addfld1d(fname='QFLX_SURF_PERROAD_XX', units='mm/s', &
         avgflag='A', long_name='URBANxx pervious road surface runoff flux', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ! --- Phase 3: Register net shortwave radiation fields ---
    ! Values are valid for the next timestep (computed during albedo calc),
    ! consistent with ELM's own FSDS family.

    ptr1d => out_netsw_roof
    call hist_addfld1d(fname='NETSW_ROOF_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx roof net shortwave radiation (next-timestep)', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_netsw_sunwall
    call hist_addfld1d(fname='NETSW_SUNWALL_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx sunlit wall net shortwave radiation (next-timestep)', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_netsw_shadwall
    call hist_addfld1d(fname='NETSW_SHADWALL_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx shaded wall net shortwave radiation (next-timestep)', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_netsw_improad
    call hist_addfld1d(fname='NETSW_IMPROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx impervious road net shortwave radiation (next-timestep)', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

    ptr1d => out_netsw_perroad
    call hist_addfld1d(fname='NETSW_PERROAD_XX', units='W/m2', &
         avgflag='A', long_name='URBANxx pervious road net shortwave radiation (next-timestep)', &
         ptr_col=ptr1d, set_nourb=0._r8, c2l_scale_type='urbanf', default='active')

  !-----------------------------------------------------------------------
  subroutine urbanxx_hist_update(num_urbanl, num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Call URBANxx getters and scatter results into the column-indexed
    ! output arrays registered during urbanxx_hist_init.
    ! No-op when urbanxx_output=.false. or num_urbanl=0.
    !
    ! Snow-layer indexing: URBANxx 2D index j=1 is the top (surface) layer;
    ! ELM history levsno index nlevsno is the top layer.
    ! Mapping: out_sno_*(c, nlevsno - j + 1) = getter_2d(l_packed, j)
    !
    use ColumnType    , only : col_pp
    use column_varcon , only : icol_roof, icol_road_imperv, icol_road_perv, &
                               icol_sunwall, icol_shadewall
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer        , intent(in) :: filter_urbanc(:)

    logical :: isLL
    integer(c_int) :: status
    integer(c_int), dimension(2) :: size2D_urb, size2D_soil, size2D_sno
    integer :: idx_roof, idx_sunwall, idx_shadwall, idx_improad, idx_perroad
    integer :: fc, c, j, l, num_active

    real(c_double), allocatable, target :: buf_urb (:)  ! flat: num_urbanl * nlevurb
    real(c_double), allocatable, target :: buf_soil(:)  ! flat: num_urbanl * nlevgrnd
    real(c_double), allocatable, target :: buf_sno (:)  ! flat: num_urbanl * nlevsno
    real(c_double), allocatable, target :: buf_1d  (:)  ! flat: num_urbanl

    ! Per-surface 2D temporaries (l_packed, j)
    real(c_double), allocatable :: roof_t(:,:), sunw_t(:,:), shadw_t(:,:)
    real(c_double), allocatable :: impr_t(:,:), perr_t(:,:)
    real(c_double), allocatable :: perr_liq(:,:), perr_ice(:,:)
    real(c_double), allocatable :: roof_sno_t(:,:), impr_sno_t(:,:), perr_sno_t(:,:)
    real(c_double), allocatable :: roof_sno_liq(:,:), impr_sno_liq(:,:), perr_sno_liq(:,:)
    real(c_double), allocatable :: roof_sno_ice(:,:), impr_sno_ice(:,:), perr_sno_ice(:,:)
    real(c_double), allocatable :: roof_sno_dz(:,:),  impr_sno_dz(:,:),  perr_sno_dz(:,:)
    real(c_double), allocatable :: roof_h2osno(:), impr_h2osno(:)
    real(c_double), allocatable :: roof_snol(:), impr_snol(:), perr_snol(:)
    real(c_double), allocatable :: t_canyon_buf(:), t_bld_buf(:)
    ! Phase 2 local buffers
    real(c_double), allocatable :: eflx_sg_roof(:), eflx_sg_sunwall(:), eflx_sg_shadwall(:)
    real(c_double), allocatable :: eflx_sg_improad(:), eflx_sg_perroad(:)
    real(c_double), allocatable :: xmf_roof(:), xmf_improad(:), xmf_perroad(:)
    real(c_double), allocatable :: netlw_roof(:), netlw_sunwall(:), netlw_shadwall(:)
    real(c_double), allocatable :: netlw_improad(:), netlw_perroad(:)
    real(c_double), allocatable :: fsh_roof(:), fsh_sunwall(:), fsh_shadwall(:)
    real(c_double), allocatable :: fsh_improad(:), fsh_perroad(:)
    real(c_double), allocatable :: qflxsnomelt_roof(:), qflxsnomelt_improad(:), qflxsnomelt_perroad(:)
    real(c_double), allocatable :: qflxsubsnow_roof(:), qflxsubsnow_improad(:), qflxsubsnow_perroad(:)
    real(c_double), allocatable :: qflxsurf_roof(:), qflxsurf_improad(:), qflxsurf_perroad(:)

    if (.not. urbanxx_output) return
    if (num_urbanl == 0) return

    isLL = UrbanKokkosIsLayoutLeft()
    size2D_urb (1) = num_urbanl;  size2D_urb (2) = nlevurb
    size2D_soil(1) = num_urbanl;  size2D_soil(2) = nlevgrnd
    size2D_sno (1) = num_urbanl;  size2D_sno (2) = nlevsno

    allocate(buf_urb (num_urbanl * nlevurb))
    allocate(buf_soil(num_urbanl * nlevgrnd))
    allocate(buf_sno (num_urbanl * nlevsno))
    allocate(buf_1d  (num_urbanl))

    allocate(roof_t(num_urbanl, nlevurb));   allocate(sunw_t(num_urbanl, nlevurb))
    allocate(shadw_t(num_urbanl, nlevurb));  allocate(impr_t(num_urbanl, nlevgrnd))
    allocate(perr_t(num_urbanl, nlevgrnd))
    allocate(perr_liq(num_urbanl, nlevgrnd)); allocate(perr_ice(num_urbanl, nlevgrnd))
    allocate(roof_sno_t(num_urbanl, nlevsno)); allocate(impr_sno_t(num_urbanl, nlevsno))
    allocate(perr_sno_t(num_urbanl, nlevsno))
    allocate(roof_sno_liq(num_urbanl, nlevsno)); allocate(impr_sno_liq(num_urbanl, nlevsno))
    allocate(perr_sno_liq(num_urbanl, nlevsno))
    allocate(roof_sno_ice(num_urbanl, nlevsno)); allocate(impr_sno_ice(num_urbanl, nlevsno))
    allocate(perr_sno_ice(num_urbanl, nlevsno))
    allocate(roof_sno_dz(num_urbanl, nlevsno));  allocate(impr_sno_dz(num_urbanl, nlevsno))
    allocate(perr_sno_dz(num_urbanl, nlevsno))
    allocate(roof_h2osno(num_urbanl)); allocate(impr_h2osno(num_urbanl))
    allocate(roof_snol(num_urbanl));   allocate(impr_snol(num_urbanl))
    allocate(perr_snol(num_urbanl))
    allocate(t_canyon_buf(num_urbanl)); allocate(t_bld_buf(num_urbanl))
    allocate(eflx_sg_roof(num_urbanl));     allocate(eflx_sg_sunwall(num_urbanl))
    allocate(eflx_sg_shadwall(num_urbanl)); allocate(eflx_sg_improad(num_urbanl))
    allocate(eflx_sg_perroad(num_urbanl))
    allocate(xmf_roof(num_urbanl));    allocate(xmf_improad(num_urbanl))
    allocate(xmf_perroad(num_urbanl))
    allocate(netlw_roof(num_urbanl));    allocate(netlw_sunwall(num_urbanl))
    allocate(netlw_shadwall(num_urbanl)); allocate(netlw_improad(num_urbanl))
    allocate(netlw_perroad(num_urbanl))
    allocate(fsh_roof(num_urbanl));    allocate(fsh_sunwall(num_urbanl))
    allocate(fsh_shadwall(num_urbanl)); allocate(fsh_improad(num_urbanl))
    allocate(fsh_perroad(num_urbanl))
    allocate(qflxsnomelt_roof(num_urbanl)); allocate(qflxsnomelt_improad(num_urbanl))
    allocate(qflxsnomelt_perroad(num_urbanl))
    allocate(qflxsubsnow_roof(num_urbanl)); allocate(qflxsubsnow_improad(num_urbanl))
    allocate(qflxsubsnow_perroad(num_urbanl))
    allocate(qflxsurf_roof(num_urbanl)); allocate(qflxsurf_improad(num_urbanl))
    allocate(qflxsurf_perroad(num_urbanl))

    ! ==============================================================
    ! Fetch all quantities from URBANxx
    ! ==============================================================

    ! --- Phase 1A: layer temperatures ---
    call UrbanGetLayerTempRoof(urbanxx, c_loc(buf_urb), size2D_urb, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_urb, roof_t, num_urbanl, nlevurb, isLL)

    call UrbanGetLayerTempSunlitWall(urbanxx, c_loc(buf_urb), size2D_urb, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_urb, sunw_t, num_urbanl, nlevurb, isLL)

    call UrbanGetLayerTempShadedWall(urbanxx, c_loc(buf_urb), size2D_urb, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_urb, shadw_t, num_urbanl, nlevurb, isLL)

    call UrbanGetLayerTempImperviousRoad(urbanxx, c_loc(buf_soil), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_soil, impr_t, num_urbanl, nlevgrnd, isLL)

    call UrbanGetLayerTempPerviousRoad(urbanxx, c_loc(buf_soil), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_soil, perr_t, num_urbanl, nlevgrnd, isLL)

    ! --- Phase 1B: snow state scalars ---
    call UrbanGetH2OSnoRoof(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    roof_h2osno(:) = buf_1d(:)

    call UrbanGetH2OSnoImperviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    impr_h2osno(:) = buf_1d(:)

    call UrbanGetNumberOfActiveSnowLayersRoof(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    roof_snol(:) = buf_1d(:)

    call UrbanGetNumberOfActiveSnowLayersImperviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    impr_snol(:) = buf_1d(:)

    call UrbanGetNumberOfActiveSnowLayersPerviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    perr_snol(:) = buf_1d(:)

    ! --- Phase 1C: snow layer fields ---
    call UrbanGetSnowLayerTempRoof(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, roof_sno_t, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerTempImperviousRoad(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, impr_sno_t, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerTempPerviousRoad(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, perr_sno_t, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerLiquidWaterRoof(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, roof_sno_liq, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerLiquidWaterImperviousRoad(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, impr_sno_liq, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerLiquidWaterPerviousRoad(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, perr_sno_liq, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerIceContentRoof(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, roof_sno_ice, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerIceContentImperviousRoad(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, impr_sno_ice, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerIceContentPerviousRoad(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, perr_sno_ice, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerThicknessRoof(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, roof_sno_dz, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerThicknessImperviousRoad(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, impr_sno_dz, num_urbanl, nlevsno, isLL)

    call UrbanGetSnowLayerThicknessPerviousRoad(urbanxx, c_loc(buf_sno), size2D_sno, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_sno, perr_sno_dz, num_urbanl, nlevsno, isLL)

    ! --- Phase 1D: pervious road soil water ---
    call UrbanGetSoilLiquidWaterPerviousRoad(urbanxx, c_loc(buf_soil), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_soil, perr_liq, num_urbanl, nlevgrnd, isLL)

    call UrbanGetSoilIceContentPerviousRoad(urbanxx, c_loc(buf_soil), size2D_soil, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    call reshape_flat_to_2d(buf_soil, perr_ice, num_urbanl, nlevgrnd, isLL)

    ! --- Phase 1E: canyon air and building temperatures ---
    call UrbanGetCanyonAirTemperature(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    t_canyon_buf(:) = buf_1d(:)

    call UrbanGetBuildingTemperature(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    t_bld_buf(:) = buf_1d(:)

    ! --- Phase 2: flux and energy balance ---
    call UrbanGetEflxSoilGrndRoof(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    eflx_sg_roof(:) = buf_1d(:)

    call UrbanGetEflxSoilGrndSunlitWall(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    eflx_sg_sunwall(:) = buf_1d(:)

    call UrbanGetEflxSoilGrndShadedWall(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    eflx_sg_shadwall(:) = buf_1d(:)

    call UrbanGetEflxSoilGrndImperviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    eflx_sg_improad(:) = buf_1d(:)

    call UrbanGetEflxSoilGrndPerviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    eflx_sg_perroad(:) = buf_1d(:)

    call UrbanGetXmfRoof(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    xmf_roof(:) = buf_1d(:)

    call UrbanGetXmfImperviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    xmf_improad(:) = buf_1d(:)

    call UrbanGetXmfPerviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    xmf_perroad(:) = buf_1d(:)

    call UrbanGetNetLongwaveRoof(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netlw_roof(:) = buf_1d(:)

    call UrbanGetNetLongwaveSunlitWall(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netlw_sunwall(:) = buf_1d(:)

    call UrbanGetNetLongwaveShadedWall(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netlw_shadwall(:) = buf_1d(:)

    call UrbanGetNetLongwaveImperviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netlw_improad(:) = buf_1d(:)

    call UrbanGetNetLongwavePerviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netlw_perroad(:) = buf_1d(:)

    call UrbanGetSensibleHeatFluxRoof(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    fsh_roof(:) = buf_1d(:)

    call UrbanGetSensibleHeatFluxSunlitWall(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    fsh_sunwall(:) = buf_1d(:)

    call UrbanGetSensibleHeatFluxShadedWall(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    fsh_shadwall(:) = buf_1d(:)

    call UrbanGetSensibleHeatFluxImperviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    fsh_improad(:) = buf_1d(:)

    call UrbanGetSensibleHeatFluxPerviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    fsh_perroad(:) = buf_1d(:)

    call UrbanGetQflxSnowMeltRoof(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    qflxsnomelt_roof(:) = buf_1d(:)

    call UrbanGetQflxSnowMeltImperviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    qflxsnomelt_improad(:) = buf_1d(:)

    call UrbanGetQflxSnowMeltPerviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    qflxsnomelt_perroad(:) = buf_1d(:)

    call UrbanGetQflxSubSnowRoof(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    qflxsubsnow_roof(:) = buf_1d(:)

    call UrbanGetQflxSubSnowImperviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    qflxsubsnow_improad(:) = buf_1d(:)

    call UrbanGetQflxSubSnowPerviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    qflxsubsnow_perroad(:) = buf_1d(:)

    call UrbanGetQflxSurfRoof(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    qflxsurf_roof(:) = buf_1d(:)

    call UrbanGetQflxSurfImperviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    qflxsurf_improad(:) = buf_1d(:)

    call UrbanGetQflxSurfPerviousRoad(urbanxx, c_loc(buf_1d), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    qflxsurf_perroad(:) = buf_1d(:)

    ! ==============================================================
    ! Scatter packed getter results into column-indexed output arrays.
    ! Non-urban columns remain at their initialized spval.
    ! ==============================================================

    idx_roof     = 0
    idx_sunwall  = 0
    idx_shadwall = 0
    idx_improad  = 0
    idx_perroad  = 0

    do fc = 1, num_urbanc
      c = filter_urbanc(fc)

      select case (col_pp%itype(c))

      case (icol_roof)
        idx_roof = idx_roof + 1
        l = idx_roof
        do j = 1, nlevurb
          out_tsoi_roof(c, j) = real(roof_t(l, j), r8)
        end do
        out_h2osno_roof(c) = real(roof_h2osno(l), r8)
        out_snol_roof(c)   = real(roof_snol(l),   r8)
        ! j=1 = top layer → ELM levsno index nlevsno
        num_active = nint(roof_snol(l))
        do j = 1, num_active
          out_sno_t_roof  (c, nlevsno - j + 1) = real(roof_sno_t  (l, j), r8)
          out_sno_liq_roof(c, nlevsno - j + 1) = real(roof_sno_liq(l, j), r8)
          out_sno_ice_roof(c, nlevsno - j + 1) = real(roof_sno_ice(l, j), r8)
          out_sno_dz_roof (c, nlevsno - j + 1) = real(roof_sno_dz (l, j), r8)
        end do
        out_t_canyon(c) = real(t_canyon_buf(l), r8)
        out_t_bld(c)    = real(t_bld_buf(l),    r8)
        out_eflx_sg_roof    (c) = real(eflx_sg_roof    (l), r8)
        out_xmf_roof        (c) = real(xmf_roof        (l), r8)
        out_netlw_roof      (c) = real(netlw_roof       (l), r8)
        out_fsh_roof        (c) = real(fsh_roof         (l), r8)
        out_qflxsnomelt_roof(c) = real(qflxsnomelt_roof (l), r8)
        out_qflxsubsnow_roof(c) = real(qflxsubsnow_roof (l), r8)
        out_qflxsurf_roof   (c) = real(qflxsurf_roof    (l), r8)

      case (icol_sunwall)
        idx_sunwall = idx_sunwall + 1
        l = idx_sunwall
        do j = 1, nlevurb
          out_tsoi_sunwall(c, j) = real(sunw_t(l, j), r8)
        end do
        out_eflx_sg_sunwall(c) = real(eflx_sg_sunwall(l), r8)
        out_netlw_sunwall  (c) = real(netlw_sunwall   (l), r8)
        out_fsh_sunwall    (c) = real(fsh_sunwall      (l), r8)

      case (icol_shadewall)
        idx_shadwall = idx_shadwall + 1
        l = idx_shadwall
        do j = 1, nlevurb
          out_tsoi_shadwall(c, j) = real(shadw_t(l, j), r8)
        end do
        out_eflx_sg_shadwall(c) = real(eflx_sg_shadwall(l), r8)
        out_netlw_shadwall  (c) = real(netlw_shadwall   (l), r8)
        out_fsh_shadwall    (c) = real(fsh_shadwall      (l), r8)

      case (icol_road_imperv)
        idx_improad = idx_improad + 1
        l = idx_improad
        do j = 1, nlevgrnd
          out_tsoi_improad(c, j) = real(impr_t(l, j), r8)
        end do
        out_h2osno_improad      (c) = real(impr_h2osno      (l), r8)
        out_snol_improad        (c) = real(impr_snol         (l), r8)
        out_eflx_sg_improad     (c) = real(eflx_sg_improad   (l), r8)
        out_xmf_improad         (c) = real(xmf_improad       (l), r8)
        out_netlw_improad       (c) = real(netlw_improad      (l), r8)
        out_fsh_improad         (c) = real(fsh_improad        (l), r8)
        out_qflxsnomelt_improad (c) = real(qflxsnomelt_improad(l), r8)
        out_qflxsubsnow_improad (c) = real(qflxsubsnow_improad(l), r8)
        out_qflxsurf_improad    (c) = real(qflxsurf_improad   (l), r8)
        num_active = nint(impr_snol(l))
        do j = 1, num_active
          out_sno_t_improad  (c, nlevsno - j + 1) = real(impr_sno_t  (l, j), r8)
          out_sno_liq_improad(c, nlevsno - j + 1) = real(impr_sno_liq(l, j), r8)
          out_sno_ice_improad(c, nlevsno - j + 1) = real(impr_sno_ice(l, j), r8)
          out_sno_dz_improad (c, nlevsno - j + 1) = real(impr_sno_dz (l, j), r8)
        end do

      case (icol_road_perv)
        idx_perroad = idx_perroad + 1
        l = idx_perroad
        do j = 1, nlevgrnd
          out_tsoi_perroad   (c, j) = real(perr_t  (l, j), r8)
          out_soilliq_perroad(c, j) = real(perr_liq(l, j), r8)
          out_soilice_perroad(c, j) = real(perr_ice(l, j), r8)
        end do
        out_snol_perroad        (c) = real(perr_snol         (l), r8)
        out_eflx_sg_perroad     (c) = real(eflx_sg_perroad   (l), r8)
        out_xmf_perroad         (c) = real(xmf_perroad        (l), r8)
        out_netlw_perroad       (c) = real(netlw_perroad       (l), r8)
        out_fsh_perroad         (c) = real(fsh_perroad         (l), r8)
        out_qflxsnomelt_perroad (c) = real(qflxsnomelt_perroad(l), r8)
        out_qflxsubsnow_perroad (c) = real(qflxsubsnow_perroad(l), r8)
        out_qflxsurf_perroad    (c) = real(qflxsurf_perroad   (l), r8)
        num_active = nint(perr_snol(l))
        do j = 1, num_active
          out_sno_t_perroad  (c, nlevsno - j + 1) = real(perr_sno_t  (l, j), r8)
          out_sno_liq_perroad(c, nlevsno - j + 1) = real(perr_sno_liq(l, j), r8)
          out_sno_ice_perroad(c, nlevsno - j + 1) = real(perr_sno_ice(l, j), r8)
          out_sno_dz_perroad (c, nlevsno - j + 1) = real(perr_sno_dz (l, j), r8)
        end do

      end select
    end do

    deallocate(buf_urb, buf_soil, buf_sno, buf_1d)
    deallocate(roof_t, sunw_t, shadw_t, impr_t, perr_t)
    deallocate(perr_liq, perr_ice)
    deallocate(roof_sno_t, impr_sno_t, perr_sno_t)
    deallocate(roof_sno_liq, impr_sno_liq, perr_sno_liq)
    deallocate(roof_sno_ice, impr_sno_ice, perr_sno_ice)
    deallocate(roof_sno_dz,  impr_sno_dz,  perr_sno_dz)
    deallocate(roof_h2osno, impr_h2osno)
    deallocate(roof_snol, impr_snol, perr_snol)
    deallocate(t_canyon_buf, t_bld_buf)
    deallocate(eflx_sg_roof, eflx_sg_sunwall, eflx_sg_shadwall, eflx_sg_improad, eflx_sg_perroad)
    deallocate(xmf_roof, xmf_improad, xmf_perroad)
    deallocate(netlw_roof, netlw_sunwall, netlw_shadwall, netlw_improad, netlw_perroad)
    deallocate(fsh_roof, fsh_sunwall, fsh_shadwall, fsh_improad, fsh_perroad)
    deallocate(qflxsnomelt_roof, qflxsnomelt_improad, qflxsnomelt_perroad)
    deallocate(qflxsubsnow_roof, qflxsubsnow_improad, qflxsubsnow_perroad)
    deallocate(qflxsurf_roof, qflxsurf_improad, qflxsurf_perroad)

  end subroutine urbanxx_hist_update

  !-----------------------------------------------------------------------
  subroutine urbanxx_hist_update_radiation(num_urbanl, num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Extract URBANxx net shortwave fields after urbanxx_netShortwave runs
    ! inside the doalb block.  Values represent the radiation computed for
    ! the *next* timestep, consistent with ELM's own shortwave history fields.
    ! No-op when urbanxx_output=.false. or num_urbanl=0.
    !
    use ColumnType    , only : col_pp
    use column_varcon , only : icol_roof, icol_road_imperv, icol_road_perv, &
                               icol_sunwall, icol_shadewall
    !
    implicit none
    integer(c_int), intent(in) :: num_urbanl
    integer(c_int), intent(in) :: num_urbanc
    integer        , intent(in) :: filter_urbanc(:)

    integer(c_int) :: status
    integer :: idx_roof, idx_sunwall, idx_shadwall, idx_improad, idx_perroad
    integer :: fc, c, l
    real(c_double), allocatable, target :: buf(:)  ! flat: num_urbanl
    real(c_double), allocatable :: netsw_roof(:), netsw_sunwall(:), netsw_shadwall(:)
    real(c_double), allocatable :: netsw_improad(:), netsw_perroad(:)

    if (.not. urbanxx_output) return
    if (num_urbanl == 0) return

    allocate(buf(num_urbanl))
    allocate(netsw_roof(num_urbanl));    allocate(netsw_sunwall(num_urbanl))
    allocate(netsw_shadwall(num_urbanl)); allocate(netsw_improad(num_urbanl))
    allocate(netsw_perroad(num_urbanl))

    call UrbanGetNetShortwaveRoof(urbanxx, c_loc(buf), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netsw_roof(:) = buf(:)

    call UrbanGetNetShortwaveSunlitWall(urbanxx, c_loc(buf), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netsw_sunwall(:) = buf(:)

    call UrbanGetNetShortwaveShadedWall(urbanxx, c_loc(buf), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netsw_shadwall(:) = buf(:)

    call UrbanGetNetShortwaveImperviousRoad(urbanxx, c_loc(buf), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netsw_improad(:) = buf(:)

    call UrbanGetNetShortwavePerviousRoad(urbanxx, c_loc(buf), num_urbanl, status)
    if (status /= URBAN_SUCCESS) call UrbanError(0, __LINE__, status)
    netsw_perroad(:) = buf(:)

    idx_roof     = 0
    idx_sunwall  = 0
    idx_shadwall = 0
    idx_improad  = 0
    idx_perroad  = 0

    do fc = 1, num_urbanc
      c = filter_urbanc(fc)
      select case (col_pp%itype(c))
      case (icol_roof)
        idx_roof = idx_roof + 1;  l = idx_roof
        out_netsw_roof(c) = real(netsw_roof(l), r8)
      case (icol_sunwall)
        idx_sunwall = idx_sunwall + 1;  l = idx_sunwall
        out_netsw_sunwall(c) = real(netsw_sunwall(l), r8)
      case (icol_shadewall)
        idx_shadwall = idx_shadwall + 1;  l = idx_shadwall
        out_netsw_shadwall(c) = real(netsw_shadwall(l), r8)
      case (icol_road_imperv)
        idx_improad = idx_improad + 1;  l = idx_improad
        out_netsw_improad(c) = real(netsw_improad(l), r8)
      case (icol_road_perv)
        idx_perroad = idx_perroad + 1;  l = idx_perroad
        out_netsw_perroad(c) = real(netsw_perroad(l), r8)
      end select
    end do

    deallocate(buf, netsw_roof, netsw_sunwall, netsw_shadwall, netsw_improad, netsw_perroad)

  end subroutine urbanxx_hist_update_radiation
    ! Reshape a Kokkos-layout flat buffer to 2D (nl, nlev).
    implicit none
    real(c_double), intent(in)  :: flat(:)
    real(c_double), intent(out) :: arr2d(:,:)
    integer        , intent(in)  :: nl, nlev
    logical        , intent(in)  :: isLayoutLeft

    integer :: l, j, count

    count = 0
    if (isLayoutLeft) then
      do j = 1, nlev
        do l = 1, nl
          count = count + 1
          arr2d(l, j) = flat(count)
        end do
      end do
    else
      do l = 1, nl
        do j = 1, nlev
          count = count + 1
          arr2d(l, j) = flat(count)
        end do
      end do
    end if

  end subroutine reshape_flat_to_2d

end module UrbanxxHistMod
