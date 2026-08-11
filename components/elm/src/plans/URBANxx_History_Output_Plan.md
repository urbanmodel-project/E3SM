# ELM–URBANxx: History Output Plan (`urbanxx_output`)

## Implementation status

| Phase | Status | Notes |
|---|---|---|
| Namelist plumbing (`urbanxx_output` flag) | ✅ Done | `elm_varctl.F90`, `controlMod.F90` |
| `UrbanxxHistMod.F90` — new module | ✅ Done | `biogeophys/UrbanxxHistMod.F90` |
| `elm_initializeMod.F90` wiring | ✅ Done | `urbanxx_hist_init` called after `elm_inst_biogeophys` |
| `elm_driver.F90` wiring | ✅ Done | `urbanxx_hist_update` called after `urbanxx_updateSnowFraction` |
| Phase 1A — layer temperatures (`levurb`, `levgrnd`) | ✅ Done | 5 fields |
| Phase 1B — snow state scalars (1D) | ✅ Done | 5 fields; `SNOWDP_ROOF_XX` deferred (requires derivation) |
| Phase 1C — snow layer state (`levsno`) | ✅ Done | 12 fields |
| Phase 1D — pervious road soil water (`levgrnd`) | ✅ Done | 2 fields |
| Phase 1E — canyon/building temperatures (1D, roof column) | ✅ Done | 2 fields |
| Phase 2 — flux and energy balance fields | ✅ Done | 27 fields; `ERRSOI_ROOF_XX` deferred (no getter) |
| Phase 3 — radiation fields (`doalb` call site) | ✅ Done | `urbanxx_hist_update_radiation` after `urbanxx_netShortwave_check` |

---


## Motivation

URBANxx parity checks log divergences to `lnd.log` but write nothing to ELM's
history NetCDF files. To characterize the known order-of-operations divergence
between ELM and URBANxx (particularly the `nstep=203` pervious-road blocker) and
to publish URBANxx results alongside ELM, we need URBANxx's own computed state
and fluxes available for post-processing.

With `urbanxx_always_continue = .true.` the run already completes without stopping.
This feature adds the companion capability: saving URBANxx data to history output.

---

## Design decisions

### Flag: `urbanxx_output` (logical, default `.false.`)

A new namelist variable gates all memory allocation, field registration, and
per-timestep extraction. When `.false.` (the default), no arrays are allocated,
no `hist_addfld` calls are made, and `urbanxx_hist_update` is a no-op — zero
overhead for production runs.

The flag is independent of `use_urbanxx`. Setting `urbanxx_output = .true.`
requires `use_urbanxx = .true.` (guarded in the init call).

### New Fortran module: `UrbanxxHistMod.F90`

Lives in `biogeophys/` alongside the other `Urbanxx*.F90` modules. Owns:
- All column-indexed output arrays (module-level allocatables)
- `urbanxx_hist_init` — allocates arrays + registers `hist_addfld` fields
- `urbanxx_hist_update` — calls URBANxx getters and fills the arrays

### Naming convention: `{ELM_BASENAME}_{SURFACE}_XX`

The `_XX` suffix marks the field as coming from URBANxx (not ELM's own solver).
Surface qualifiers: `ROOF`, `SUNWALL`, `SHADWALL`, `IMPROAD`, `PERROAD`.

Examples: `TSOI_ROOF_XX`, `H2OSNO_IMPROAD_XX`, `SNO_T_PERROAD_XX`.

Fields that are landunit-level in URBANxx (canyon air temperature, building
temperature) have no surface qualifier and are assigned to the `icol_roof`
column (one roof column per urban landunit, providing a clean 1-to-1 mapping):
`T_CANYON_XX`, `T_BLD_XX`.

### Memory layout

Output arrays are allocated `(begc:endc)` for 1D fields and
`(begc:endc, 1:nlev)` for 2D fields, matching the shape `hist_addfld1d`/
`hist_addfld2d` expect from `ptr_col`. Non-urban columns are initialized to
`spval` and never written; urban-specific fields use `set_nourb=0._r8` and
`c2l_scale_type='urbanf'` (same as ELM's own urban fields like `URBAN_AC`).

Estimated memory per global run at ~100K columns (worst case):
- 1D field: ~800 KB
- `levurb` (5-layer) 2D field: ~4 MB
- `levgrnd` (15-layer) 2D field: ~12 MB
- `levsno` (5-layer) 2D field: ~4 MB
- ~40 fields total: roughly 100–150 MB when enabled, 0 when not

### l-to-c scatter pattern

URBANxx getters return flat `(num_urbanl)` or `(num_urbanl * nlev)` arrays,
indexed by packed landunit position. The scatter back to ELM's column-indexed
output arrays reuses the same pattern as the existing parity checks: a single
loop over `filter_urbanc` with per-surface-type incrementing counters
(`idx_roof`, `idx_sunwall`, `idx_shadwall`, `idx_improad`, `idx_perroad`), using
`col_pp%itype(c)` to route each column to the right packed index.

```fortran
idx_roof = 0; idx_sunwall = 0; idx_shadwall = 0
idx_improad = 0; idx_perroad = 0
do fc = 1, num_urbanc
  c = filter_urbanc(fc)
  select case (col_pp%itype(c))
  case (icol_roof)
    idx_roof = idx_roof + 1
    out_h2osno_roof(c) = getter_h2osno_roof(idx_roof)
    do j = 1, nlevurb
      out_tsoi_roof(c, j) = getter_tsoi_roof_2d(idx_roof, j)
    end do
    ! canyon air temperature and building temperature assigned here too
    out_t_canyon(c)  = getter_t_canyon(idx_roof)
    out_t_bld(c)     = getter_t_bld(idx_roof)
  case (icol_sunwall)
    idx_sunwall = idx_sunwall + 1
    do j = 1, nlevurb
      out_tsoi_sunwall(c, j) = getter_tsoi_sunwall_2d(idx_sunwall, j)
    end do
  ! ... icol_shadewall, icol_road_imperv, icol_road_perv ...
  end select
end do
```

The 2D reshape from flat to `(num_urbanl, nlev)` follows the existing pattern
in `UrbanxxSoilTemperatureMod.F90`:
```fortran
do l = 1, num_urbanl
  do j = 1, nlevurb
    getter_tsoi_roof_2d(l, j) = getter_buf_roof((l-1)*nlevurb + j)
  end do
end do
```

### Extraction call placement

A single `urbanxx_hist_update` call after `urbanxx_updateSnowFraction`
(the last URBANxx call in the main physics block, currently line 1052 in
`elm_driver.F90`). At this point all thermodynamics, snow physics, hydrology,
and fluxes have been computed for the current timestep.

Net shortwave fields are computed inside the `if (doalb)` block (line ~1483,
runs only on radiation timesteps). Two options:
- **Option A (simpler, deferred)**: Exclude shortwave from the first
  implementation; add a second `urbanxx_hist_update_radiation` call inside
  `doalb` in Phase 2.
- **Option B**: Accept that `NETSW_{SURF}_XX` is the value computed for the
  *next* timestep (consistent with how ELM's own `FSDS` family works) and
  extract it from the `doalb` block.

**Recommendation**: Option A for Phase 1. The primary diagnostic value is in
state (temperatures, snow layers) and energy-balance fluxes, not shortwave.

---

## Field catalogue

### Phase 1: State fields (highest diagnostic value) — ✅ Done

#### 1A — Layer temperatures (2D)

| Field name | Getter | `type2d` | Notes |
|---|---|---|---|
| `TSOI_ROOF_XX` | `UrbanGetLayerTempRoof` | `levurb` | 5-layer |
| `TSOI_SUNWALL_XX` | `UrbanGetLayerTempSunlitWall` | `levurb` | |
| `TSOI_SHADWALL_XX` | `UrbanGetLayerTempShadedWall` | `levurb` | |
| `TSOI_IMPROAD_XX` | `UrbanGetLayerTempImperviousRoad` | `levgrnd` | 15-layer |
| `TSOI_PERROAD_XX` | `UrbanGetLayerTempPerviousRoad` | `levgrnd` | |

#### 1B — Snow state scalars (1D)

| Field name | Getter | Units | Notes |
|---|---|---|---|
| `H2OSNO_ROOF_XX` | `UrbanGetH2OSnoRoof` | kg/m² | SWE |
| `H2OSNO_IMPROAD_XX` | `UrbanGetH2OSnoImperviousRoad` | kg/m² | |
| `SNOWDP_ROOF_XX` | *(derived from snow layers once N>0; thin-snow scalar otherwise)* | m | ⬜ deferred |
| `SNOL_ROOF_XX` | `UrbanGetNumberOfActiveSnowLayersRoof` | count | int→real |
| `SNOL_IMPROAD_XX` | `UrbanGetNumberOfActiveSnowLayersImperviousRoad` | count | |
| `SNOL_PERROAD_XX` | `UrbanGetNumberOfActiveSnowLayersPerviousRoad` | count | |

#### 1C — Snow layer state (2D with `levsno`)

| Field name | Getter | Notes |
|---|---|---|
| `SNO_T_ROOF_XX` | `UrbanGetSnowLayerTempRoof` | Layer temp (K) |
| `SNO_T_IMPROAD_XX` | `UrbanGetSnowLayerTempImperviousRoad` | |
| `SNO_T_PERROAD_XX` | `UrbanGetSnowLayerTempPerviousRoad` | |
| `SNO_LIQ_ROOF_XX` | `UrbanGetSnowLayerLiquidWaterRoof` | kg/m² |
| `SNO_LIQ_IMPROAD_XX` | `UrbanGetSnowLayerLiquidWaterImperviousRoad` | |
| `SNO_LIQ_PERROAD_XX` | `UrbanGetSnowLayerLiquidWaterPerviousRoad` | |
| `SNO_ICE_ROOF_XX` | `UrbanGetSnowLayerIceContentRoof` | kg/m² |
| `SNO_ICE_IMPROAD_XX` | `UrbanGetSnowLayerIceContentImperviousRoad` | |
| `SNO_ICE_PERROAD_XX` | `UrbanGetSnowLayerIceContentPerviousRoad` | |
| `SNO_DZ_ROOF_XX` | `UrbanGetSnowLayerThicknessRoof` | m |
| `SNO_DZ_IMPROAD_XX` | `UrbanGetSnowLayerThicknessImperviousRoad` | |
| `SNO_DZ_PERROAD_XX` | `UrbanGetSnowLayerThicknessPerviousRoad` | |

#### 1D — Pervious road soil water (2D with `levgrnd`)

| Field name | Getter | Units |
|---|---|---|
| `SOILLIQ_PERROAD_XX` | `UrbanGetSoilLiquidWaterPerviousRoad` | kg/m² |
| `SOILICE_PERROAD_XX` | `UrbanGetSoilIceContentPerviousRoad` | kg/m² |

#### 1E — Building / canyon scalars (1D, assigned to roof column)

| Field name | Getter | Units |
|---|---|---|
| `T_CANYON_XX` | `UrbanGetCanyonAirTemperature` | K |
| `T_BLD_XX` | `UrbanGetBuildingTemperature` | K |

---

### Phase 2: Flux and energy balance fields (1D) — ✅ Done

| Field name | Getter | Units | Surfaces |
|---|---|---|---|
| `EFLX_SG_ROOF_XX` | `UrbanGetEflxSoilGrndRoof` | W/m² | ROOF |
| `EFLX_SG_SUNWALL_XX` | `UrbanGetEflxSoilGrndSunlitWall` | W/m² | SUNWALL |
| `EFLX_SG_SHADWALL_XX` | `UrbanGetEflxSoilGrndShadedWall` | W/m² | SHADWALL |
| `EFLX_SG_IMPROAD_XX` | `UrbanGetEflxSoilGrndImperviousRoad` | W/m² | IMPROAD |
| `EFLX_SG_PERROAD_XX` | `UrbanGetEflxSoilGrndPerviousRoad` | W/m² | PERROAD |
| `XMF_ROOF_XX` | `UrbanGetXmfRoof` | W/m² | ROOF |
| `XMF_IMPROAD_XX` | `UrbanGetXmfImperviousRoad` | W/m² | IMPROAD |
| `XMF_PERROAD_XX` | `UrbanGetXmfPerviousRoad` | W/m² | PERROAD |
| `ERRSOI_ROOF_XX` | derived in `UrbanxxEnergyBalanceMod` | W/m² | ROOF (after Phase 2) | ⬜ deferred — no getter |
| `NETLW_ROOF_XX` | `UrbanGetNetLongwaveRoof` | W/m² | all surfaces |
| `NETLW_SUNWALL_XX` | `UrbanGetNetLongwaveSunlitWall` | W/m² | |
| `NETLW_SHADWALL_XX` | `UrbanGetNetLongwaveShadedWall` | W/m² | |
| `NETLW_IMPROAD_XX` | `UrbanGetNetLongwaveImperviousRoad` | W/m² | |
| `NETLW_PERROAD_XX` | `UrbanGetNetLongwavePerviousRoad` | W/m² | |
| `FSH_ROOF_XX` | `UrbanGetSensibleHeatFluxRoof` | W/m² | ROOF |
| `FSH_IMPROAD_XX` | `UrbanGetSensibleHeatFluxImperviousRoad` | W/m² | |
| `FSH_PERROAD_XX` | `UrbanGetSensibleHeatFluxPerviousRoad` | W/m² | |
| `FSH_SUNWALL_XX` | `UrbanGetSensibleHeatFluxSunlitWall` | W/m² | |
| `FSH_SHADWALL_XX` | `UrbanGetSensibleHeatFluxShadedWall` | W/m² | |
| `QFLX_SNOMELT_ROOF_XX` | `UrbanGetQflxSnowMeltRoof` | mm/s | ROOF |
| `QFLX_SNOMELT_IMPROAD_XX` | `UrbanGetQflxSnowMeltImperviousRoad` | mm/s | |
| `QFLX_SNOMELT_PERROAD_XX` | `UrbanGetQflxSnowMeltPerviousRoad` | mm/s | |
| `QFLX_SUBSNOW_ROOF_XX` | `UrbanGetQflxSubSnowRoof` | mm/s | ROOF |
| `QFLX_SUBSNOW_IMPROAD_XX` | `UrbanGetQflxSubSnowImperviousRoad` | mm/s | |
| `QFLX_SUBSNOW_PERROAD_XX` | `UrbanGetQflxSubSnowPerviousRoad` | mm/s | |
| `QFLX_SURF_ROOF_XX` | `UrbanGetQflxSurfRoof` | mm/s | ROOF |
| `QFLX_SURF_IMPROAD_XX` | `UrbanGetQflxSurfImperviousRoad` | mm/s | |
| `QFLX_SURF_PERROAD_XX` | `UrbanGetQflxSurfPerviousRoad` | mm/s | |

---

### Phase 3: Radiation fields (requires `doalb` call site) — ✅ Done

| Field name | Getter | Units |
|---|---|---|
| `NETSW_ROOF_XX` | `UrbanGetNetShortwaveRoof` | W/m² |
| `NETSW_SUNWALL_XX` | `UrbanGetNetShortwaveSunlitWall` | W/m² |
| `NETSW_SHADWALL_XX` | `UrbanGetNetShortwaveShadedWall` | W/m² |
| `NETSW_IMPROAD_XX` | `UrbanGetNetShortwaveImperviousRoad` | W/m² |
| `NETSW_PERROAD_XX` | `UrbanGetNetShortwavePerviousRoad` | W/m² |

These are extracted inside `if (doalb)` in `elm_driver.F90`, from a separate
`urbanxx_hist_update_radiation` call or as an optional argument to
`urbanxx_hist_update`. The value is valid for the *next* timestep (computed
during albedo calculation), consistent with ELM's own shortwave history fields.

---

## Files to change

| File | Change |
|---|---|
| `main/elm_varctl.F90` | Add `logical, public :: urbanxx_output = .false.` next to the existing URBANxx flags |
| `main/controlMod.F90` | Add `urbanxx_output` to `elm_inparm` namelist declaration + `mpi_bcast` |
| `bld/namelist_files/namelist_defaults.xml` | `<urbanxx_output>.false.</urbanxx_output>` |
| `bld/namelist_files/namelist_definition.xml` | Description entry in the URBANxx block |
| `biogeophys/UrbanxxHistMod.F90` | **New file** — arrays, `urbanxx_hist_init`, `urbanxx_hist_update` |
| `main/elm_initializeMod.F90` | `use UrbanxxHistMod`; call `urbanxx_hist_init` inside the `use_urbanxx` init block |
| `main/elm_driver.F90` | `use UrbanxxHistMod`; call `urbanxx_hist_update` after `urbanxx_updateSnowFraction` |

No `CMakeLists.txt` change is required — CIME's build system auto-discovers all
`.F90` files in `biogeophys/`.

---

## `UrbanxxHistMod.F90` internal structure

```
module UrbanxxHistMod
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use elm_varctl   , only : urbanxx_output
  use spmdMod      , only : masterproc
  use elm_varpar   , only : nlevurb, nlevgrnd, nlevsno
  use histFileMod  , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use decompMod    , only : bounds_type
  use urban_mod    ! URBANxx C API

  implicit none
  private

  public :: urbanxx_hist_init
  public :: urbanxx_hist_update

  ! ---- module-level output arrays (allocated only when urbanxx_output) ----
  ! 2D state
  real(r8), allocatable :: out_tsoi_roof(:,:)      ! (begc:endc, nlevurb)
  real(r8), allocatable :: out_tsoi_sunwall(:,:)
  real(r8), allocatable :: out_tsoi_shadwall(:,:)
  real(r8), allocatable :: out_tsoi_improad(:,:)   ! (begc:endc, nlevgrnd)
  real(r8), allocatable :: out_tsoi_perroad(:,:)
  ! snow layer state
  real(r8), allocatable :: out_sno_t_roof(:,:)     ! (begc:endc, nlevsno)
  real(r8), allocatable :: out_sno_t_improad(:,:)
  real(r8), allocatable :: out_sno_t_perroad(:,:)
  real(r8), allocatable :: out_sno_liq_roof(:,:)
  ...
  ! 1D scalars
  real(r8), allocatable :: out_h2osno_roof(:)      ! (begc:endc)
  real(r8), allocatable :: out_snol_roof(:)
  real(r8), allocatable :: out_t_canyon(:)
  real(r8), allocatable :: out_t_bld(:)
  ...

  ! saved decomposition bounds (needed by urbanxx_hist_update)
  integer :: saved_begc, saved_endc

contains

  subroutine urbanxx_hist_init(num_urbanl, bounds)
    ! Allocates arrays and registers hist_addfld fields.
    ! Called once from elm_initializeMod under use_urbanxx .and. urbanxx_output.
  end subroutine

  subroutine urbanxx_hist_update(num_urbanl, num_urbanc, filter_urbanc)
    ! Calls URBANxx getters → scatters into column-indexed output arrays.
    ! No-op when .not. urbanxx_output.
  end subroutine

end module UrbanxxHistMod
```

### `urbanxx_hist_init` logic
1. Save `bounds%begc` / `bounds%endc` to module variables.
2. Allocate all output arrays; initialize to `spval`.
3. For each field, call `hist_addfld1d` or `hist_addfld2d` with:
   - `fname` = field name (e.g. `'TSOI_ROOF_XX'`)
   - `avgflag='A'` (time-average; use `'I'` for instantaneous layer counts)
   - `ptr_col` pointing to the module-level array
   - `set_nourb=0._r8` (zero for non-urban columns in gridcell averages)
   - `c2l_scale_type='urbanf'` (standard urban scaling)
   - `default='inactive'` for diagnostics that are rarely needed (snow layer
     geometry fields `SNO_DZ`, `SNO_LIQ`, soil water 2D); `default='active'`
     for primary state fields (temperatures, `H2OSNO`, `SNOL`).
   - Snow layer 2D fields need `no_snow_behavior=no_snow_normal`
     (same as ELM's `SNO_T`).

### `urbanxx_hist_update` logic
1. Return immediately if `.not. urbanxx_output`.
2. Allocate flat getter buffers (same `num_urbanl * nlev` sizing as existing
   modules; reuse local `allocate/deallocate` pattern).
3. Call each URBANxx getter into its flat buffer.
4. Reshape flat buffers to 2D `(num_urbanl, nlev)`.
5. Single pass over `filter_urbanc` with per-surface-type packed index counters
   to scatter into `out_*` arrays. Non-urban columns remain `spval`.
6. Deallocate flat getter buffers.

---

## Implementation order (within Phase 1)

1. **Namelist plumbing** (5 files, ~20 lines total) — enables the flag to be
   set before any code runs; validates the wiring.
2. **`UrbanxxHistMod.F90` skeleton** — allocates 1D arrays, registers them with
   `hist_addfld1d`, wires init/update calls. Start with just `H2OSNO_ROOF_XX`
   and `T_CANYON_XX` end-to-end to confirm the output appears in history NetCDF.
3. **Add remaining 1D scalars** (Phase 1B–1E, ~15 fields).
4. **Add 2D layer temperature fields** (Phase 1A, `levurb` + `levgrnd`).
5. **Add 2D snow layer fields** (Phase 1C, `levsno`); note ELM's
   `no_snow_behavior=no_snow_normal` is required so non-snow timesteps write
   `spval` rather than 0.
6. **Validate**: run `all_urban.1x1_brazil` (warm, no snow) — all snow-related
   fields should be `spval`; temperature and scalar flux fields should match ELM
   values. Then run `modified_f19_f19` with `urbanxx_always_continue = .true.`
   for a full month to collect the divergence data for post-processing.
7. **Phase 2 and 3** (flux fields and radiation) in separate commits.

---

## Known constraints and gotchas

- **`begc` must be preserved across init/update**: `hist_addfld` stores a raw
  Fortran pointer to the module-level array; the array must remain allocated
  and un-moved for the entire run. Module-level `allocatable` with saved bounds
  satisfies this.

- **`no_snow_behavior`**: Required for all `levsno` 2D fields. Without it, ELM's
  history machinery may fill non-active snow layers with 0 instead of `spval`,
  producing misleading averages.

- **`type2d='levurb'`**: Already defined in `histFileMod.F90` (line 1902). Valid
  for 5-layer urban walls and roof.

- **Snow-layer indexing convention**: URBANxx's snow layers are indexed
  `0..nlevsno-1` (most recently accumulated at index 0, deepest at
  `nlevsno-1`). ELM's `levsno` dimension in history is indexed `1..nlevsno`
  with layer 1 being the top. The scatter loop must apply a `+1` offset or
  reverse the order to match ELM's convention — confirm against ELM's `SNO_T`
  output before finalizing.

- **`errsoi` for walls**: Walls carry no water/ice state (Xmf=0 hardcoded).
  The errsoi residual for walls is available from `UrbanxxEnergyBalanceMod`'s
  local computation but is not currently exposed via a getter. Either add a
  getter in Phase 2 or compute it in-place in `urbanxx_hist_update` using
  the existing formula from `UrbanxxEnergyBalanceMod.F90`.

- **`doalb` shortwave**: `urbanxx_netShortwave` only runs when `doalb=.true.`
  (typically every N timesteps). If Phase 3 shortwave fields are registered with
  `hist_addfld`, the arrays must be updated only inside `doalb` and left
  unchanged on non-radiation timesteps — the history accumulator will then
  average the stale (previous-radiation-step) value on those steps, which is the
  same behavior as ELM's own `FSDS`. Document this clearly in the field's
  `long_name`.
