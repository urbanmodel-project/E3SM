---
name: port-elm-module
description: "Use when porting a new ELM Fortran physics module to URBANxx Kokkos C++. Provides a step-by-step checklist covering C++ implementation, C API, Fortran bindings, ELM integration, and correctness verification."
---

# Port an ELM Physics Module to URBANxx

Use this prompt when adding a new physics module to URBANxx. Follow each step in order — later steps depend on earlier ones.

## Context

- URBANxx source: `external_models/urbanxx/src/`
- C API header: `external_models/urbanxx/include/Urban.h`
- Fortran bindings: `external_models/urbanxx/include/urban_mod.F90`
- ELM integration wrappers: `biogeophys/Urbanxx*.F90`
- Central ELM driver: `biogeophys/UrbanxxMod.F90`
- ELM ↔ URBANxx index mapping: see `CLAUDE.md` (0-based C++ ↔ 1-based Fortran; separate Kokkos views per surface type)

---

## Step 1 — Implement physics in C++ with Kokkos

File: `external_models/urbanxx/src/UrbanXxx.cpp`

- Use `Kokkos::parallel_for` with a `RangePolicy` over landunit index `l = 0..nUrbanLandunits-1`
- Access per-surface data via `urban.surfaceName` (e.g., `urban.perviousRoad`, `urban.roof`)
- Avoid host-side conditionals inside kernels; use Kokkos view accessors
- Match the ELM Fortran logic exactly — cross-reference the source Fortran file while writing

## Step 2 — Expose via the C API

File: `external_models/urbanxx/include/Urban.h`

- Add setter/getter declarations following the existing naming pattern:
  - Setters: `UrbanSetXxxYyy(UrbanType* urban, int l, double val)`
  - Getters: `UrbanGetXxxYyy(UrbanType* urban, int l, double* val)`
- Add implementations in `external_models/urbanxx/src/UrbanParamsType.cpp` (params) or the relevant `Xxx.cpp`
- Add private implementation headers under `include/private/` if needed

## Step 3 — Bind to Fortran

File: `external_models/urbanxx/include/urban_mod.F90`

- Add `bind(C)` interface declarations for each new C function
- Follow the existing pattern: `use iso_c_binding`, `type(c_ptr), value :: urban`, `integer(c_int), value :: l`, etc.
- Fortran index `l_fortran` → C index: pass `int(l_fortran - 1, c_int)` (adjust for 0-based URBANxx)

## Step 4 — Write the ELM wrapper module

File: `biogeophys/UrbanxxXxxMod.F90` (new file, following `UrbanxxSoilWaterMod.F90` as template)

Key responsibilities:
- Loop over urban landunits: iterate over ELM columns/pfts that map to each URBANxx `l`
- Translate ELM column index `c` → URBANxx landunit index `l` (see `CLAUDE.md` mapping table)
- Call C API setters to push input data from ELM into URBANxx
- Call the URBANxx physics routine
- Call C API getters to pull output data back into ELM arrays
- Handle any unit conversions (ELM uses SI; verify units match)

## Step 5 — Register in the ELM driver

File: `biogeophys/UrbanxxMod.F90`

- Add `use UrbanxxXxxMod, only: urbanxx_xxx`
- Add `call urbanxx_xxx(...)` in the correct position in the timestep sequence:
  ```
  urbanxx_SetAtmosphericForcing → urbanxx_netLongwave → urbanxx_surfaceFluxes →
  urbanxx_soilTemperature → urbanxx_soilWater → urbanxx_netShortwave
  [after ELM Drainage:] urbanxx_pervRoad_drainage
  ```
- Pass the required ELM data structures (bounds, col, lun, atm2lnd_inst, etc.)

## Step 6 — Add a correctness check

Reference: `urbanxx_soilWater_check` (in `UrbanxxSoilWaterMod.F90`) as a template.

- After calling the URBANxx routine, loop over layers/surfaces and compare ELM vs URBANxx values
- Default tolerance: 1×10⁻⁸ (match existing checks)
- On mismatch: print a per-layer/per-surface table then call `endrun`
- Print format: include ELM value, URBANxx value, absolute difference, and layer/surface index

## Step 7 — Build and verify

```bash
# 1. Rebuild URBANxx
cd /Users/gautam.bisht/projects/urban/e3sm/components/elm/src/external_models/urbanxx/for_e3sm
cmake --build . --target urban -j$(sysctl -n hw.logicalcpu)

# 2. Rebuild ELM
cd /Users/gautam.bisht/projects/urban/e3sm/cime/scripts/all_urban.1x1_brazil.I1850ELM.PNNL-L07D666226.gnu11.54cd255e6e.2026-02-17
./case.build --skip-provenance-check

# 3. Run
./case.submit
```

- Check `run/lnd.log.<LID>` for any mismatch abort messages
- If the check fires: add temporary debug prints for the matrix coefficients / input arrays to locate the discrepancy (see existing `[ELM pre-Drainage]` / `[URBANxx pre-Drainage]` pattern in `SoilHydrologyMod.F90` and `UrbanHydrology.cpp`)

## Checklist

- [ ] C++ Kokkos kernel implemented and matches ELM Fortran logic
- [ ] C API setter/getter declared in `Urban.h` and implemented
- [ ] `bind(C)` Fortran interface added to `urban_mod.F90`
- [ ] ELM wrapper `UrbanxxXxxMod.F90` created
- [ ] Wrapper registered and called in `UrbanxxMod.F90` in the correct sequence position
- [ ] Correctness check added with `endrun` on failure
- [ ] Both libraries rebuilt; `case.submit` completes without mismatch abort
- [ ] Any temporary debug prints removed (or documented in `CLAUDE.md` if bug remains open)
