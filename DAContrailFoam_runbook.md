# DAFoam contrailFoam integration runbook

## Objective

Implement the smallest correct first DAFoam integration for the reduced CASSANDRA multicomponent gas-phase solver, with aerosol and PBE physics disabled.

This path intentionally targets the no-PBE configuration:

- aerosolActive false;
- massCouplingActive false;
- heatCouplingActive false;

The goal is to support the gas-phase transient compressible multicomponent solver while deferring the full aerosol/PBE and dynamic-LES turbulence work.

## Working assumptions

- CASSANDRA standalone solver: contrailFoam
- DAFoam fork: feature/contrailfoam-primal-wrapper
- OpenFOAM: v2506
- DAFoam baseline already validated in original, ADR, and ADF modes
- The reduced solver must preserve the compressible PIMPLE structure from DARhoPimpleFoam
- The species transport pattern comes from the standalone contrailFoam equations
- dynamicKEqn in the DAFoam turbulence adjoint layer is not yet supported and remains out of scope for this phase

## Scope boundaries

This implementation intentionally does not include:

- aerosol moment transport
- PBE / microphysics / source coupling
- latent vapor sink or latent heat source terms
- dynamicKEqn turbulence adjoint suppor
- incomplete DAContrailFoam stub registration without DAStateInfo and DAResidual initialization
- modifications to generated lnInclude or Make/linux* directories

## What is logically needed for this first correct implementation

The first implementation should be built around the minimum set of solver and model features needed to reproduce the standalone no-PBE gas-phase dynamics in DAFoam.

### 1. A proper compressible PIMPLE solver family

- A new DAFoam runtime-selectable solver family is required because the existing compressible solvers do not encode the CASSANDRA multicomponent gas-phase structure.
- The solver must follow the transient compressible PIMPLE workflow used by `DARhoPimpleFoam` rather than copying a constant-property simplification.
- The real OpenFOAM references are the compressible pressure-based solver family (`rhoPimpleFoam` / `rhoSimpleFoam`) and the CASSANDRA `contrailFoam` implementation.

### 2. Thermodynamic model compatibility

- The reduced solver must work with `psiReactionThermo` and `multiComponentMixture`.
- The gas model is not a chemistry-disabled reacting-flow model in the usual sense; it is a multicomponent ideal-gas transport solver with `hePsiThermo` and no net chemical source.
- The required thermodynamic closure is the OpenFOAM mixture-thermo machinery, not ad hoc constant-property reconstruction.

### 3. Species transport states and reconstruction

- Every active non-inert species must be represented as a discrete transported state in the DAFoam residual structure.
- The inert species (`N2`) must be reconstructed from the mass-fraction sum constraint, not treated as a separate transported field in the first pass.
- The residual pattern should match the standalone `YEqn_aero.H` form:
  - transient term
  - conservative convection term
  - diffusive term using `muEff/Sc`
  - zero source term in the reduced no-PBE path

### 4. Enthalpy and pressure-correction coupling

- The sensible-enthalpy equation must include the same pressure-work and kinetic-energy terms used in the standalone solver.
- The pressure correction and face-flux update must occur through the same PIMPLE structure used in the compressible pressure-based OpenFOAM solver.
- The first pass should not include aerosol-derived latent source terms or microphysics coupling.

### 5. Turbulence-model scope discipline

- DAFoam has no current support for `dynamicKEqn` in the turbulence adjoint layer.
- For the first correct implementation, the turbulence model used by the primal must be one already represented by DAFoam, or the case must be restricted to a supported model.
- This is a scope guard, not a physics limitation of the standalone solver.

### 6. Derivative and validation discipline

- The initial implementation should focus on primal consistency and residual structure.
- Only after the reduced gas-phase state is converged should the adjoint derivative path be validated.
- In the future, the full PBE/aerosol extension will require explicit treatment of positivity and source splitting, because those operations are non-smooth and must be handled carefully for adjoint work.

### 7. Reduced species block boundary

The reduced implementation should solve the active non-inert species equations in the same conservative finite-volume form used by the standalone solver,

- transient term: $\partial_t(\rho Y_k)$
- convection term: $\nabla \cdot (\phi Y_k)$
- diffusion term: $\nabla \cdot (\mu_{\mathrm{eff}}/Sc \nabla Y_k)$
- source term: zero in the no-PBE, no-aerosol stage

The inert species is then reconstructed from the mass-fraction constraint,

$$
Y_{\mathrm{inert}} = 1 - \sum_{k \in \mathcal{A}} Y_k,
$$

where $\mathcal{A}$ is the set of active non-inert species. This is the correct reduced gas-phase boundary for the first DAFoam patch; the latent-vapor and microphysics source terms remain disabled until the full aerosol/PBE model is implemented.

## What DAContrailFoam can do (current capabilities)

DAContrailFoam is a DAFoam runtime-selectable solver family
(`solverName: "DAContrailFoam"`, registered in the `Compressible` registry) that
targets the reduced CASSANDRA contrailFoam gas-phase physics — compressible,
transient, multicomponent, no aerosol/PBE. Verified and structural capabilities
after Stage 1+2:

- **Instantiable end to end in original mode:** `PYDAFOAM` construction with
  `solverName: "DAContrailFoam"` completes: the C++ solver registers, reads the
  case, and initializes every layer (verified by the container smoke test).
- **Full thermodynamics:** `psiReactionThermo` (`hePsiThermo` +
  `multiComponentMixture`), Sutherland transport, ideal-gas closure, inert
  specie `N2` validated at startup; owns and registers `p`, `T`, `thermo:psi`,
  and species fields `O2`, `CO2`, `H2O`, `N2` from the case `0/` directory.
- **Complete primal field set:** `U`, `rho`, `phi`, `dpdt`, `K`, pressure
  control, `MRFProperties` registry entry — i.e. everything standalone
  contrailFoam's `createFields.H` builds, in DAFoam ownership conventions.
- **Turbulence-agnostic selection:** reads `simulationType` from
  `constant/turbulenceProperties` and instantiates either LES (`LESModel`, the
  case's `dynamicKEqn`) or RAS — no RAS-only hard-code — plus the DAFoam
  `DATurbulenceModel` adjoint-side wrapper of the same model.
- **Full adjoint-side object graph:** `createAdjoint.H` runs and its every
  registered-object lookup succeeds (DAModel, DAStateInfo, DAIndex, DAField,
  DAResidualContrailFoam, DACheckMesh, DALinearEqn, objective function list).
- **Defined adjoint state vector:** `U, p, T, O2, CO2, H2O, nut` (per cell) and
  `phi` (per face) — count verified exact against the mesh
  (`9*nCells + nFaces = 8,337,125` on the NASA PAL case).
- **Physically validated transient primal (Stage 3, Phases A–C):**
  `solvePrimal()` runs the full outer time loop with CASSANDRA ordering
  (`readTimeControls` → `storePrevIter` → `CourantNo`/`setDeltaT` → `++runTime`
  → `rhoEqnContrail` → PIMPLE(`UEqn`→`YEqn`→`hEqn`→`pEqn`) → turbulence
  `correct()` → `rho = thermo.rho()` → write), with `reduceIO` intermediate
  writes via `writeAdjStates` and a full `runTime.write()` on the last step.
  `DADynamicKEqn` k-residual dimensions fixed in Phase A so
  `calcPrimalResidualStatistics` passes (`Total Residual Norm2 12447191.93`).
  Field-level validation (Phase C) against standalone contrailFoam: step-1
  gates bit-exact except the documented FP-level k source seed; long-run
  snapshot at t=0.06 agrees on means to ≤3.5% (`UMean`) / ≤1.0% (`TMean`,
  `rhoMean`), inlet patch values bit-identical; outlet-BC investigation found
  no boundary-treatment defect (see Phase C findings below).

In short: discovery, construction, field creation, turbulence selection,
adjoint plumbing, the physics time loop, and primal field validation are
complete in original mode. The remaining work is the Phase D adjoint path:
time-averaged mean-observation objective → ADR enablement audit → ADR/ADF
rebuild (explicitly gated) → finite-difference-validated derivatives.

## Progress so far

### Step 1: solver family skeleton

Created the first minimal DAFoam solver family skeleton for a dedicated DAContrailFoam entry:

- src/adjoint/DASolver/DAContrailFoam/DAContrailFoam.H
- src/adjoint/DASolver/DAContrailFoam/DAContrailFoam.C

This follows the lifecycle of DARhoPimpleFoam and keeps the solver runtime-selectable but intentionally simple.

### Step 2: registration

Updated the DAFoam solver registry and build list so the new solver is discoverable:

- src/adjoint/Make/files
- dafoam/pyDAFoam.py

This is the minimal build-level requirement for the new runtime-selectable solver family.

### Step 3: DAStateInfo registration

Created a dedicated DAStateInfo class:

- src/adjoint/DAStateInfo/DAStateInfoContrailFoam.H
- src/adjoint/DAStateInfo/DAStateInfoContrailFoam.C

It registers the reduced gas-phase state family and preserves the solver metadata expected by DAFoam:

- p
- T
- Y (multicomponent species aggregate state)
- U
- phi
- nut

The key design decision was to do this without pretending to support the full aerosol/PBE system yet.

### Step 4: DAResidual skeleton

Created the corresponding DAResidual class:

- src/adjoint/DAResidual/DAResidualContrailFoam.H
- src/adjoint/DAResidual/DAResidualContrailFoam.C

This class follows the same compressible transient structure as the rhoPimple residual, but the multicomponent species transport is kept to the reduced no-PBE form.

### Step 5: reduced multicomponent species residual path

Extended the residual implementation to include the reduced no-PBE species-transport logic based on the standalone contrailFoam formulation:

- solve active non-inert species mass fractions
- reconstruct the inert species as $1 - \sum Y_k$
- use the effective diffusion form with $\mu_{eff}/Sc$
- omit aerosol/PBE source terms

This matches the physical intent of the reduced solver while still staying within a minimal first implementation.

## Current status

Phases A–C of the validation plan are complete (details in the sections below):

- **Phase A — residual validation: passed.** kResidual dimensions fixed to
  `rhoDimensions()*k_.dimensions()/dimTime` (DADynamicKEqn.C); smoke test and
  `--residuals` both pass.
- **Phase B / Stage 3 — primal parity: passed.** Time loop + equation parity
  committed as `6f7a7ca` "Add Stage 3 primal time loop and contrail equation
  parity" (9 files, +288/−104, incl. `rhoEqnContrail.H`).
- **Phase C — standalone-vs-DAFoam comparison: executed; engineering-level
  pass, formally partial.** Gates + 0.06 s twin runs + figure export done;
  the 0.3 s run was stopped early (t≈0.076) for the outlet-BC investigation,
  which found no BC defect (below). The remaining optional action is the
  end-of-run (0.3 s) and time-averaged-field comparison if a full-length run
  is wanted; Phase D has been authorized to start in parallel.
- **Phase D — D1 complete, D2 implemented and validated (uncommitted).** D1
  design in `DAContrailFoam_phaseD1_mean_objective_design.md`; D2 delivered the
  `functionTimeSteps_` bounds guard (`DASolver.C` `calcAllFunctions` store),
  the `--func` driver option that injects `meanTOutlet`, an original-mode
  rebuild with a clean log, and a 0.012 s validation run against `fieldAverage1`
  (`evalFunctions` 208.0155615012 K vs `TMean` outlet areaAverage
  208.0159243 K — difference −1.74e-6 rel, entirely dt-weighting). Full write-up
  in the **Phase D2** section below; session state in
  `DAContrailFoam_session_handover.md`.

## Recommended next step

Phase D per `dafoam_contrailfoam_copilot_prompt_pack.md` (state in
`DAContrailFoam_session_handover.md`, resume prompts in its §3):

1. **D2 — done (see the Phase D2 section).** Bounds guard added and rebuilt
   (clean log); `meanTOutlet` injected via `stage3_primal_compare.py --func`;
   0.012 s validation vs `fieldAverage1` completed and the dt-weighting
   difference quantified (−1.74e-6 rel). Remaining D2 follow-ups are optional:
   a dt-weighted `timeOp` if exact agreement with `fieldAverage` is wanted,
   and recompute-from-t0 on restart before any adjoint use.
2. **D3** — one-control/one-objective sensitivity + identifiability plan
   (no ADR).
3. **D4** — ADR enablement audit (checklist + build/test plan); only then,
   on explicit request, rebuild ADR/ADF and validate derivatives by finite
   difference on one short-window mean objective.

## Files to watch

- cassandra/applications/solvers/contrailFoam/createFields.H
- cassandra/applications/solvers/contrailFoam/YEqn_aero.H
- cassandra/applications/solvers/contrailFoam/hEqn_aero.H
- dafoam/src/adjoint/DAResidual/DAResidualRhoPimpleFoam.C
- dafoam/src/adjoint/DAStateInfo/DAStateInfoRhoPimpleFoam.C
- dafoam/src/adjoint/DASolver/DARhoPimpleFoam/DARhoPimpleFoam.C

## Step 6 (Stage 1): original-mode compile fixes

Two root compile errors were fixed in the DAContrailFoam translation unit (the rest
of the file's diagnostics were cascade warnings from these):

1. `createRefsContrail.H` — `fv::convectionScheme<scalar>::New(...).ptr()` called
   `tmp<>::ptr()`, which is deleted in OF v2506 (the "clone" cascade in `tmpI.H`).
   Fixed to keep the returned `tmp` directly:

   ```cpp
   tmp<fv::convectionScheme<scalar>> mvConvection =
       fv::convectionScheme<scalar>::New(mesh, fields, phi, mesh.divScheme("div(phi,Yi_h)"));
   ```

   This is the same expression `DAResidualContrailFoam.C` (compiles) and CASSANDRA
   `YEqn_aero.H` use; consumers use `mvConvection->fvmDiv(...)`, which `tmp` supports.

2. `YEqnContrail.H` — `volScalarField sumY("sumY", mesh, dimensionedScalar(...))`
   used a constructor that no longer exists in v2506. Fixed to the CASSANDRA
   `YEqn_aero.H` pattern: `volScalarField sumY(0.0*Y[0]);`.

## Step 7 (Stage 2): initSolver / createFields structural fix

The skeleton's `createFieldsContrail.H` looked up a `psiReactionThermo` that was
never created and re-created `p`/`T`/`phi` in the wrong registry instance. It was
rewritten so every object the residual/adjoint machinery looks up is created and
registered, in the correct order, matching standalone CASSANDRA field creation and
DAFoam ownership conventions (`DARhoPimpleFoam`):

- `createFieldsContrail.H` (rewritten), executed inside `initSolver()`:
  1. `pThermoPtr_.reset(psiReactionThermo::New(mesh))` — creates/registers
     `thermophysicalProperties` and owns `p`, `T`, `thermo:psi`, and the species
     `Y` fields (all `MUST_READ` from the start-time directory);
     `thermo.validate(args.executable(), "h", "e")`;
     inert-specie validation copied from CASSANDRA `createFields.H`
     (`FatalIOErrorIn` if `inertSpecie` not in the species list);
  2. `rho` with `READ_IF_PRESENT`/`AUTO_WRITE` initialized from `thermo.rho()`;
     `U` `MUST_READ`; `phi` via `compressibleCreatePhiPython.H`
     (declares `phi` and sets `phiPtr_`, DAFoam's standard include);
     `mesh.setFluxRequired(p.name())`;
  3. `pressureControlPtr_` constructed with `pimple.dict()` (proven by standalone
     CASSANDRA and DAPimpleFoam);
  4. `compressible::turbulenceModel::New(rho, U, phi, thermo)` — registers `nut`
     (an adjoint state) and `alphat` (a residual lookup), both present in the
     generated `0/` directory;
  5. `dpdt` (dimensions `p/dimTime`) and `K = 0.5*magSqr(U)` — the two derived
     fields `DAResidualContrailFoam` looks up, same as `DARhoPimpleFoam`;
  6. `MRFPtr_.reset(new IOMRFZoneListDF(mesh))` — registers `MRFProperties`
     (reads the dict only if present; the case has none).
- `DAContrailFoam.H` — `pThermoPtr_` changed from raw pointer to
  `autoPtr<psiReactionThermo>` (thermo owns `p`/`T`, so the old `pPtr_`/`TPtr_`
  members were removed to avoid double ownership/double registration); added
  `pressureControlPtr_` and `MRFPtr_`; added `#include "pressureControl.H"`.
- `DAContrailFoam.C::initSolver()` — reordered to the proven DAPimpleFoam shape:
  `args` → `createPimpleControlPython.H` → `createFieldsContrail.H` →
  turbulence model name read from `constant/turbulenceProperties` with an
  **LES-or-RAS branch** (`LES` subdict `LESModel` if present, else `RASModel`;
  no RAS-only hard-code, per the prompt-pack stop condition) →
  `DATurbulenceModel::New(...)` (must precede `DAModel`) → `#include "createAdjoint.H"`
  (previously missing entirely; creates `DAModel`, `DAStateInfo`, `DAIndex`,
  `DAField`, `DAResidual`, `DACheckMesh`, `DALinearEqn`, function list).
- `createRefsContrail.H` — reordered so `thermo` is bound first;
  `volScalarField& p = thermo.p()` (not a dangling `pPtr_`); the unused, null
  `T` dereference removed (`T` only ever used as `thermo.T()`); everything else
  (`rho`, `U`, `phi`, `K`, `dpdt`, `psi`, species, `MRF`, `turbulence`) now
  dereferences pointers created in Step 7.

Consistency argument for each choice: field creation list and inert-specie check
copy standalone CASSANDRA `createFields.H`; pointer ownership, `validate`,
`thermo.rho()`, `dpdt`/`K` construction, pressure control and the initSolver
sequence copy `DARhoPimpleFoam`/`DAPimpleFoam`; equation files were left byte-
identical to their CASSANDRA ports (`UEqnContrail`/`YEqnContrail`/`hEqnContrail`/
`pEqnContrail`).

## Verification performed

1. **Compile (original mode):** detached build via
   `scripts/build_dafoam_source.sh`; log `logs/dafoam-Allmake-20260929-152848.log`
   shows `Ctoo: DAContrailFoam.C`, `link: libDASolver.so`, `Build Successful!`,
   **0 errors**, one expected warning (`unused runTime` in `createRefs`, kept for
   Stage 3's time loop). The log, not `dafoam-Allmake-status.txt`, is the
   authoritative check (the status file can report stale SUCCEEDED).
2. **Runtime smoke test (in the `dafoam-dev` container, not on the host):**
   generated the case `0/` with `generate_inlet_BCs.py`, built the mesh with
   `blockMesh` from the case's existing `system/blockMeshDict`, then ran
   `scripts/smoke_test_dacontrail.py`. The script bypasses ADR/ADF initialization
   (the AD libs are stale: `libDASolverADR.so` predates this patch,
   `libDASolverADF.so` has no `DAContrailFoam`) and runs original mode only.
   Result **PASS**: full `PYDAFOAM` construction (= `initSolver` →
   `createFields` → turbulence/`DATurbulenceModel` → `createAdjoint`, whose every
   registered-object lookup would fatal if a field were missing) succeeded;
   adjoint state count exactly `9*nCells + nFaces = 8,337,125`; `dpdt`, `K`,
   `rho` readable through the public API (`rho=0.2708` consistent with thermo).

## Known remaining items

- **Phase C formal closure (optional):** the 0.3 s twin runs (`da_long`,
  `of_long`) were stopped at t≈0.076 for the outlet-BC investigation; a full
  run + final-field/time-mean comparison can be re-run later if a formal
  Phase C sign-off beyond the current evidence is wanted.
- **`reduceIO` intermediate writes omit `N2`, `nut`, `alphat`, `rho`** (and
  the `uniform/time` state dict); final-step writes are full. Comparison
  scripts must not expect those fields in intermediate time dirs.
- **LES `nut`/`k` pointwise divergence is real and expected:** the dynamic
  Germano model amplifies the FP-level step-1 k source seed
  (`maxRelDiff(nut)=1.88` already at step 3); comparisons need engineering
  tolerances on means, not rtol=1e-6 on instantaneous pointwise values.
- `DAField.C:1137` hard-codes `turbDict.subDict("RAS")`, but only fires when
  `primalBC` contains `useWallFunction`; default `primalBC = {}` keeps it dormant.
- ADR/ADF builds must be rebuilt before reverse/forward AD runs can use
  `DAContrailFoam` (original-mode `COMPILE_DAFOAM_NOAD=1` skips them);
  rebuild gated on the Phase D4 audit + explicit request.
- **Uncommitted in working tree:** runbook/prompt-pack status updates, the D1
  design doc (`DAContrailFoam_phaseD1_mean_objective_design.md`), and the
  session handover (`DAContrailFoam_session_handover.md`) from this session
  (checkpoint `6f7a7ca` itself is committed).

## Step 8 (Phase A): kResidual dimensions fix

`DADynamicKEqn.C` k-equation residual used `dimVolume*rhoDimensions()` where
the k-equation terms carry `rhoDimensions()*k_.dimensions()/dimTime`
(`[1 -1 -3]`, not `[1 2 -3]`). Fixed the `kRes` dimensions accordingly; smoke
test and `--residuals` now pass (`Total Residual Norm2 12447191.93`).

## Step 9 (Stage 3 / Phase B): primal time loop and equation parity

Committed as `6f7a7ca` (9 files, +288/−104). `solvePrimal` gained the outer
time loop with CASSANDRA ordering, `rhoEqnContrail.H`, ddtCorr fix in
`pEqnContrail.H`, and full write handling. Parity notes verified against the
reference case:

- `fvSolution/PIMPLE`: `consistent no`, `nOuterCorrectors 1`,
  `nCorrectors 2` — so standalone's `pimple.consistent()` branch
  (`pcEqn.H`) is never taken; DA's always-`pEqnContrail.H` is equivalent.
- `pEqnContrail.H` ports the non-transonic branch of CASSANDRA `pEqn.H`;
  transonic/pcEqn branches are rejected by design (matches
  `DAResidualContrailFoam`, documented in the file header).
- `rho = thermo.rho()` runs after the PIMPLE loop in both solvers.
- OF writes with `runTime.write()` every step (Time::write checks
  `writeTime`); DA uses `writeAdjStates` (reduceIO) for intermediate steps
  and `runTime.write()` on the last step.

## Phase C validation findings

### Gate runs (t = 1e-5 and 3e-5)

- Step 1: `p T O2 CO2 H2O rho U nut` bit-exact; `k` maxRelDiff 3.45e-6
  (energy-weighted `k/TKEsum` 2.23e-9). Source of the k seed: matrix bisect
  with instrumented kEqn dumps showed **diag bit-identical (0 diff)**, all
  dumped inputs bit-identical (k_int, rho_int, flux_int, inlet-boundary
  ν_eff/ρ/nut/k); only the explicit source vector differs (148,393 cells,
  max 1.56e-8, worst cell 210) — FP-level, not a structural bug.
- Step 3: `p T rho` species OK; `U` 5.75e-6, `k` 1.35e-2, `nut` 1.88 — the
  dynamic Germano model amplifies the seed (expected chaotic LES behavior).

### Twin 0.06 s runs (`da_long` vs `of_long`, endTime 0.3, writeInterval 0.06)

Both started from identical `0/` (byte-identical `diff -rq`), ran to
t≈0.076 and were stopped for the investigation; both wrote the t=0.06 step.
Comparison at t=0.06 (`scripts/bc_stats_006.py`):

| region | U | p | T | k |
|---|---|---|---|---|
| inlet patch | **0 (bit-exact)** | — | **0 (bit-exact)** | **0 (bit-exact)** |
| inlet-adjacent cells (rel) | 2.0e-5 | 2.0e-5 | 1.3e-5 | 4.6e-5 |
| outlet-adjacent cells (rel) | 0.168 | 2.3e-3 | 0.042 | 0.121 |
| worst cells anywhere (rel) | 0.306 | 2.5e-3 | 0.099 | 0.455 |

- `phi` flux balance: inlet sum matches to 1.9e-5 rel
  (−0.0813250 vs −0.0813265); outlet sum da 0.088664 vs of 0.085166 (+4.1%);
  outlet patch faces differ up to 79% face-wise but are consistently signed
  (da: faster/hotter/lower-p exhaust).
- Time-averaged fields (0→0.06, `fieldAverage`): `UMean` rel 3.5%,
  `TMean` 0.8%, `rhoMean` 1.0%; second moments (`UPrime2Mean` etc.) 5–10% —
  consistent with a short window + chaotic divergence.

### Outlet boundary-condition investigation (user-reported suspicion)

Verdict: **no outlet-BC treatment defect.** Evidence:

1. `0/` directories byte-identical (all BC types/values shared);
   `system/controlDict` differs only in `startFrom`; `constant/
   turbulenceProperties` differs only by DA's required `turbulence on;`.
2. Loop/call-order parity with standalone (equation order, `rho = thermo.rho()`
   placement, `correctBoundaryConditions` via the same OF field machinery);
   inlet patch values at t=0.06 are bit-identical — the solvers apply BCs
   through identical code and identical dictionaries.
3. Divergence grows monotonically from the inlet (FP-level) to the outlet;
   the worst differences are in the plume interior, not at the outlet. All
   outlet biases are mutually physically consistent (hotter + faster +
   lower-p exhaust = different plume arrival), i.e. the outlet BC
   (pressureInletOutletVelocity + waveTransmissive) is **transmitting** a
   diverged interior state, not generating it.
4. At step 1 the outlet patch values are identical by the same argument that
   gave bit-exact internal fields (identical adjacent cells + identical
   patch dicts + identical phi/dt evaluated by the same BC classes).

### VTK output naming (`da_long_1` vs `of_long_1685` — cosmetic)

Both directories are t=0.06 (`.vtm.series` `"time": 0.06` for both).
`foamToVTK` names output dirs `_ + runTime.timeIndex()` (foamToVTK.C:722),
and `Time::setTime` overwrites `timeIndex` from the per-time-state dict
`<time>/uniform/time` → `index` (Time.C:936) when present:

- standalone OF full write → `0.06/uniform/time` with `index 1685`
  (= solver iteration count at t=0.06) → `of_long_1685`;
- DA reduceIO intermediate write (`writeAdjStates`, per-field `.write()`)
  omits `uniform/time` → fallback to foamToVTK's directory enumeration →
  `da_long_1`.

So one suffix is the saved solver iteration (OF) and the other is the
converted-directory counter (DA) — exactly a write-metadata difference, not
a controlDict or physics difference.

## Phase D2 — `meanTOutlet` mean objective (implemented + validated)

### Modified files (in this changeset)

| file | change |
|---|---|
| `src/adjoint/DASolver/DASolver.C` (`calcAllFunctions`, store formerly at :369) | bounds guard: if `listIndex >= functionTimeSteps_[idxI].size()` grow to `listIndex+1` with zero-fill before the store; skip the store entirely when `listIndex < 0` |
| `work/scripts/stage3_primal_compare.py` (host `/Users/lbosc/dafoam-source-docker/work/scripts/`, not a git repo) | new `--func` flag (mode `da` only, rejects other modes): injects `daOptions["function"]["meanTOutlet"]` = `{type: patchMean, source: patchToFace, patches: [outlet], varName: T, varType: scalar, index: 0, scale: 1.0, timeOp: average, nStepsFrac: 1.0}` and prints `evalFunctions` after `solvePrimal` as `STAGE3: func meanTOutlet = <%.12e>` |

Rebuild (original mode, `COMPILE_DAFOAM_NOAD=1`): `DASolver.C` recompiled and
`libDASolver.so` relinked; authoritative log
`/workspace/logs/dafoam-Allmake-20261001-062527.log` ends with
`*** Build Successful! ***` and contains no `error` lines.

### Window / reset / finalization / restart behavior (as built)

- **Window.** `timeOp average` → `DASolver::getTimeOpRange` (`DASolver.C:442`):
  trailing window of `max(2, round(nStepsFrac*(listEnd+1)))` **steps** ending
  at `listEnd`. `nStepsFrac = 1.0` ⇒ `startIdx = 0` ⇒ the whole recorded
  history `[0, primalFinalTimeIndex_-1]`. An absolute `[t0, t1]` window is
  emulated with `endTime = t1`, `nStepsFrac = (t1-t0)/endTime` (D1 §8).
  Caveat: the window is step-count based, not time based — with
  `adjustTimeStep` (validation run: dt 1.2e-5 → 5.15e-5) a step fraction is
  not the same fraction of elapsed time; exact time windows need fixed dt or a
  time-based/dt-weighted range.
- **Reset.** `functionTimeSteps_` is allocated exactly once, in
  `setDAFunctionList()` (called from `createAdjoint.H`, i.e. at solver
  construction), sized `round(endTime/deltaT)` and zero-filled. There is no
  window-level reset and no serialization to disk. A second `solvePrimal` in
  the same process overwrites indices `0..newFinal-1`; `getTimeOpFuncVal` only
  ever reads up to `primalFinalTimeIndex_-1`, so stale tail entries are never
  observed.
- **Finalization.** `primalFinalTimeIndex_` is set at the end of
  `DAContrailFoam::solvePrimal` (`DAContrailFoam.C:208`) and on the
  residual/funcStd early-exit in `DASolver::loop` (`DASolver.C:226`).
  `getTimeOpFuncVal` is a pure re-read of the stored window → idempotent, no
  value to freeze; steps past the final index do not exist, so nothing can
  modify the reported mean after the run.
- **Restart.** The accumulator is in-memory only. On a restart with
  `startTime > 0` the new guard keeps the store in bounds (timeIndex starts at
  the restart index), but entries `0..startIndex-1` remain zero, so a
  `nStepsFrac = 1.0` mean would be **diluted by pre-restart zeros**.
  Recompute-from-t0 / replay (D1 §3) is still not implemented and is required
  before any adjoint (D3/D4) use.
- **Bounds guard rationale.** The list is sized from the *initial* deltaT;
  a dt that shrinks after start, or a nonzero restart timeIndex, overruns it
  (unguarded `scalarList::operator[]` write → heap corruption). The guard
  grows on demand with zero-fill; `listIndex < 0` (timeIndex 0) is skipped.

### Validation run (short copied case, `fieldAverage1` active)

Case `/tmp/stage3/da_func` = copy of `da_long/{0,constant,system}` with
`endTime 0.012`, `writeInterval 0.012` (so the final step is a writeTime);
`adjustTimeStep yes`, `maxCo 0.4`; `fieldAverage1` untouched (`timeStart 0`,
`executeInterval 1`, `writeControl writeTime`, `base time`).
Driver: `stage3_primal_compare.py --mode da --case /tmp/stage3/da_func --out
/tmp/stage3/da_func.npz --func`. Result: 266 steps, `primalFail=0`, 230 s
execution, no errors, `0.012/TMean` written.

- **`TMean` boundary values confirmed present**: `boundaryField` on all three
  patches, `type calculated`, `value nonuniform List<scalar> 2800` →
  fieldAverage *does* average boundary fields, so a patch comparison is valid.
- Ground truth via `postProcess -func tmeanOutlet -time 0.012` with
  `system/tmeanOutlet` = `{type surfaceFieldValue, regionType patch,
  names (outlet), operation areaAverage, fields (TMean), writeFields false,
  log true}` → `areaAverage(outlet) of TMean = 208.0159243 K`.

| quantity | value (K) |
|---|---|
| `evalFunctions` `meanTOutlet` — unweighted step mean (`DATimeOpAverage`) | **208.0155615012** |
| fieldAverage `TMean` outlet areaAverage — dt-weighted | **208.0159243** |
| difference (DA − OF) | **−3.628e-4 K = −1.744e-6 rel (−0.0002 %)** |

**The gap is entirely dt-weighting.** Both means recomputed offline from the
per-step `meanTOutlet:` and `Time =` lines of the same run:

- unweighted mean of the 266 per-step values = `208.0155614996` → agrees with
  `evalFunctions` to 7.6e-12 rel (Python→C++ plumbing and the `nStepsFrac=1.0`
  window are exact);
- dt-weighted mean (`dt_n = t_n − t_{n-1}`, min 1.2e-5 / max 5.15e-5) =
  `208.0159242540` → agrees with fieldAverage to 2.2e-10 rel, i.e. **identical
  sample set and identical weights**: both sample the post-step states
  `t1..tN` exactly once (no `t=0` sample — `funcObj.start()` does not
  average), fieldAverage weighting each sample by `dt_n`;
- unweighted − dt-weighted from the same per-step values = `−3.6275e-4 K`,
  matching the observed `evalFunctions − fieldAverage` gap to 5e-8 K.

So the only systematic discrepancy vs `fieldAverage` is the missing dt-weighting
in `DATimeOpAverage`: **1.74e-6 relative (3.63e-4 K) for this 0–0.012 s
window**. It is ~3 orders of magnitude inside the D1 acceptance gate (Phase C
`TMean` twin spread 0.8 %). If exact bit-level agreement with fieldAverage is
wanted, add a dt-weighted `timeOp` (the one-line option flagged in D1 §8);
`evalFunctions` itself is consistent with the unweighted definition it
implements.

### Guard-trigger run (`/tmp/stage3/da_guard`)

To exercise the resize path the case was set up so the pre-sized list is
smaller than the executed step count: `endTime 0.0006` + `maxDeltaT 5e-6` ⇒
list sized `round(0.0006/1e-5) = 60`, but **120 steps executed**, so the guard
ran on steps 61–120. The run completed (`primalFail=0`, no errors), the final
`evalFunctions` = `208.0000000000` equals the offline mean of all 120 stored
values, and every printed running average — including post-resize steps — read
back exactly, i.e. no out-of-bounds write or read occurred. Caveat: the outlet
`T` is still the uniform 208 K initial state that early, so the observable is
constant in this window; the run proves survival + correct readback across the
resize, not value diversity. The main `da_func` run (266 steps vs 1200
pre-sized) did not need the guard.

## Phase D2 hardening — Prompt 2 error handling (committed)

### What fails loudly now (previously silently wrong)

| case | pre-fix behavior | post-fix |
|---|---|---|
| `getTimeOpFuncVal` with `primalFinalTimeIndex_ == 0` (primal never ran / primal failed) | window `[0, -1]` → `avg /= 0` → **NaN returned silently** | `FatalError` "primalFinalTimeIndex_ == 0 … run solvePrimal first" |
| `getdFScaling` with `primalFinalTimeIndex_ == 0` | window check never matches → **0 returned silently** (wrong scaling) | same `FatalError` in `getdFScaling` |
| `getTimeOpFuncVal("wrongName")` | loop no-match → **0.0 returned silently** (`getdFScaling` was already fatal) | `FatalError` listing the configured function names |
| `getdFScaling("wrongName", …)` | already fatal (message did not include the name) | fatal, now includes the offending name |
| empty window `iEnd < iStart` inside `DATimeOpAverage::compute` | `avg /= 0` → NaN (and `dFScaling` → `1/0` = inf) | `FatalError` with iStart/iEnd/list size |
| empty window in `DATimeOpFinal::compute` | `valList[-1]` **out-of-bounds read** | same `FatalError` guard |
| empty window in `DATimeOpMax::compute` | KS mode `log(0)` = −inf, orig mode bogus `-1e16` | same `FatalError` guard |
| `DAFunctionPatchMean` with zero total face area | `value / areaSum_` → **NaN** | `FatalError` "Total face area … is 0" |

The empty-window guards in the three `DATimeOp` classes are backstops: every
externally reachable path (a `getTimeOpFuncVal`/`getdFScaling` call with an
empty window) is already blocked by the two `primalFinalTimeIndex_ == 0`
fatals, and `calcAllFunctions` only ever runs with `timeIndex >= 1`
(`listIndex >= 0`), for which `getTimeOpRange` always returns
`startIdx <= endIdx`. The guards make the invariant explicit for future
callers.

### Files changed

`src/adjoint/DASolver/DASolver.C` (`getTimeOpFuncVal`, `getdFScaling`),
`src/adjoint/DATimeOp/DATimeOpAverage.C` (both methods),
`src/adjoint/DATimeOp/DATimeOpFinal.C` (`compute`),
`src/adjoint/DATimeOp/DATimeOpMax.C` (`compute`),
`src/adjoint/DAFunction/DAFunctionPatchMean.C` (`calcFunction`).

Rebuild: `logs/dafoam-Allmake-20261001-073832.log` — all five files
recompiled (`Ctoo: …`), `libDASolver.so` relinked (07:38), `*** Build
Successful! ***`, `grep -c error` = 0.

### Verification (container `/tmp/stage3`, case `da_err` = `da_func` copy with
`endTime 0.00002`)

Script `prompt2_test.py <mode> <case>`; each failure case is its own process
because `abort(FatalError)` kills it (exit 134 = SIGABRT).

| mode | action | exit | observed |
|---|---|---|---|
| `noprim_val` | init with `function`, **no** primal → `evalFunctions({})` | 134 | `FOAM FATAL ERROR: primalFinalTimeIndex_ == 0 … getTimeOpFuncVal` |
| `noprim_scale` | init, no primal → `solver.getdFScaling("meanTOutlet", 0)` | 134 | same message from `getdFScaling` |
| `bogus_val` | 2-step primal → `solver.getTimeOpFuncVal("bogusFunc")` | 134 | `functionName "bogusFunc" not found in daFunctionPtrList_. Configured functions: ["meanTOutlet"]` |
| `bogus_scale` | 2-step primal → `solver.getdFScaling("bogusFunc", 0)` | 134 | `functionName "bogusFunc" not found in daFunctionPtrList_.` |
| `areasum` | function dict with `source: allCells` (leaves `faceSources_` empty) → fatal inside the primal loop | 134 | `Total face area of patchMean function "meanTOutlet" is 0 …` |
| sanity | normal path: `stage3_primal_compare.py --mode da --case da_err --func` | 0 | `primalFail=0`, `meanTOutlet = 2.080000000000e+02` (2 pristine steps, expected) |

The empty-window `DATimeOp*` guards could not be triggered from Python (the
DASolver-level fatals fire first), so they are verified by code-path analysis
above rather than by a run.

### New finding (follow-up, not fixed here)

An empty patch list (`"patches": []` with `source: patchToFace`) segfaults
during `PYDAFOAM` **construction** — PETSc reports `signal 11 SEGV`, then
`MPI_ABORT` (exit 59) — before `initSolver()` completes and before any of the
guards above are reachable (none of them run during init; only the unchanged
`DAFunctionPatchMean` constructor runs there). `-X faulthandler` and
`PETSC_OPTIONS=-no_signal_trap` did not yield a Python backtrace and no
`gdb`/`catchsegv` exists in the container, so the exact frame is unisolated.
Recorded as a follow-up; the `areaSum_ == 0` case was therefore verified via
`source: allCells` instead, which reaches `calcFunction` normally.

## Notes

The implementation should stay intentionally narrow and reviewable. The first milestone is not a full DAFoam contrail solver; it is the reduced gas-phase path that matches the validated CASSANDRA no-PBE physics.
