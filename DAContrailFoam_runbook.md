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
- **Phase D — D1/D2/D3/D4 documents complete; D2 + Prompt 2 code committed.**
  Commits: `f33eab6` (D2 bounds guard + Phase A–D2 record), `fef5cc5`
  (Prompt 2 fail-fast error handling, verified with five failure-case runs),
  `dc28c07` (status docs). D1 design in
  `DAContrailFoam_phaseD1_mean_objective_design.md`; D2 delivered the
  `functionTimeSteps_` bounds guard (`DASolver.C` `calcAllFunctions` store),
  the `--func` driver option that injects `meanTOutlet`, an original-mode
  rebuild with a clean log, and a 0.012 s validation run against `fieldAverage1`
  (`evalFunctions` 208.0155615012 K vs `TMean` outlet areaAverage
  208.0159243 K — difference −1.74e-6 rel, entirely dt-weighting). D3 plan in
  `DAContrailFoam_phaseD3_sensitivity_plan.md` (inlet jet velocity scale α,
  fixed-dt FD protocol, R1–R4 checks, G1–G7 acceptance, costed run matrix);
  D4 audit in `DAContrailFoam_phaseD4_adr_audit.md` (non-smooth inventory,
  coded-BC blocker, replay/time-average requirements, G1–G6 build/test gates).
  Session state in `DAContrailFoam_session_handover.md`.

## Recommended next step

Phase D per `dafoam_contrailfoam_copilot_prompt_pack.md` (state in
`DAContrailFoam_session_handover.md`, resume prompts in its §3):

1. **D2 + Prompt 2 — done.** Bounds guard + fail-fast error handling
   (`fef5cc5`), rebuilt clean, verified by five failure-case runs (see the
   Phase D2 / Prompt 2 sections); `meanTOutlet` validated vs `fieldAverage1`
   (−1.74e-6 rel, dt-weighting only). Optional D2 follow-ups: a dt-weighted
   `timeOp` if exact `fieldAverage` agreement is wanted; recompute-from-t0 on
   restart before any adjoint use.
2. **D3 — plan done** (`DAContrailFoam_phaseD3_sensitivity_plan.md`).
   Execution is the next *compute* step: pilot run (fixed dt = 1e-5, Co
   check), R1 determinism gate, then the Tier-1 FD matrix (≈ 8.6 h serial),
   report dJ/dα with CI under G1–G7. No solver code required (control = one
   constant in `0/U`).
3. **D4 — audit done** (`DAContrailFoam_phaseD4_adr_audit.md`). The
   ADR/ADF rebuild remains gated on explicit request; when requested, start
   at gate G1 (ADR compile + coded-BC runtime compile) and finish at G5
   (FD validation of the 0.012 s `meanTOutlet` window) before any "ADR
   support" claim.

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

## Phase D3 Tier-1 execution — run results (2026-10-01)

Plan: `DAContrailFoam_phaseD3_sensitivity_plan.md` §5–§8. Executed in
`dafoam-dev` under `/tmp/stage3/d3_*`; analysis script
`scripts/d3_sensitivity_analysis.py` (`--dir`, `--t0`, `--t1`).

### Configuration and deviations from the plan

| parameter | plan §5 | executed | reason |
|---|---|---|---|
| mesh | 20/30/100 = 280,000 cells | **14/20/60 = 78,960 cells** (`generate_blockMeshDict.py`) | user-approved coarsening; 3.5× faster, identical for all runs (FD self-consistency unaffected) |
| time step | fixed `deltaT 1e-5` | fixed **`deltaT 1.25e-5`** | Co-limited by the xy cell (3.86 mm): `max\|U\|`=97.07 → Co_max = **0.3146 ≤ 0.4**; 1.25e-5 divides both window boundaries exactly |
| end time / window | 0.06 s, [0.02, 0.06] | unchanged | — |
| `nStepsFrac` | 0.6667 | 0.6667 → **last 3200 of 4800 samples** (t = 0.0200125…0.06) | plan's 4000/6000 counts assumed dt = 1e-5; window bounds unchanged |
| run set | J0 ×2, ±1 %, ±2 % | 6 matrix runs + `d3_j0s` (series cross-check) + `d3_t2p` probe to 0.15 s (diagnostic) | — |

Rate ≈ 0.6 s/step with 6 runs contended (≈ 17–30 min/run vs 86 min planned);
7 runs × ~600 MB fit the 8 GB Docker VM; no thermal warnings on the host
throughout (battery ≈ 30.9 °C, load ≈ 10/10 during the runs).

### J values (window [0.02, 0.06] s, `evalFunctions`, `primalFail = 0` everywhere)

| case | α (m/s) | J (K) |
|---|---|---|
| `d3_j0a` | 96.4 | 208.1350321309 |
| `d3_j0b` | 96.4 | 208.1350321309 (identical) |
| `d3_j0s` | 96.4 | 208.1350321309 (identical; + FO) |
| `d3_p01` | 97.364 (+1 %) | 208.1372536936 |
| `d3_m01` | 95.436 (−1 %) | 208.1322759207 |
| `d3_p02` | 98.328 (+2 %) | 208.1409050374 |
| `d3_m02` | 94.472 (−2 %) | 208.1296164240 |

### Gates (plan §8)

| gate | result | evidence |
|---|---|---|
| G1 determinism | **PASS** | `j0a` vs `j0b` final fields bit-identical (`--mode cmp`, maxAbsDiff = 0 on all fields); window series `j0a == j0b == j0s` array-equal; J strings equal |
| G2 dt sequence | **PASS** | identical `Time =` series, 4800 steps, all 7 runs |
| G3 plumbing | **PASS** | offline window mean vs `evalFunctions`: rel **5.704e-13** (≤1e-9); sample count 3200 = T_win/dt exactly; `nStepsFrac 0.6667` matches window |
| G4 stationarity | **FAIL** | half-means 208.0661450 / 208.2038084, diff **1.376e-01** vs σ_J = 4.834e-02 (calibrated z = 1.42) |
| G5 SNR | **FAIL** | \|ΔJ\| = 4.978e-03 (h=0.01) and 1.129e-02 (h=0.02) vs threshold 5√2σ_J = **3.418e-01** |
| G6 step-size | **FAIL** (formally, via G5) | raw slopes +0.2489 / +0.2822 per α̂ — same sign, dev 13.4 % (<25 %), but neither h is G5-valid |
| G7 primal health | **PASS** | `primalFail=0` ×7, series finite, Co_max = 0.3146 ≤ 0.4 |

Window statistics (baseline, [0.02, 0.06], N = 3200): σ_inst = 9.083e-02 K,
τ_int = 5.67 ms, σ_J = 4.834e-02 K.

FD table (plan window, with 95 % CI from plan §7 — CI dominated by the
non-stationarity-inflated σ_J, hence not interpretable):

| h | dJ/d(α̂) | dJ/dα (K per m/s) | curvature C |
|---|---|---|---|
| 0.01 | +0.2489 ± 6.701 | +0.002582 ± 0.0695 | −5.346 K |
| 0.02 | +0.2822 ± 3.351 | +0.002928 ± 0.0348 | +1.143 K |

### Diagnosis (why G4–G6 failed)

5 ms bin means of the baseline series: flat **plateau 208.026 K for
t = 0.005–0.030** (acoustic equilibration only), then the plume-arrival ramp
208.026 → 208.34 through the end of the run — the run stops *during* arrival,
so the window is dominated by drift, which inflates σ_inst/τ_int/σ_J ~100×.

Diagnostic probe `d3_t2p` (identical baseline, endTime 0.15, 12000 steps,
`primalFail=0`, J = 208.3224566202 over [0, 0.15]): arrival peaks near
t ≈ 0.065, then the outlet mean oscillates with **period ≈ 0.07 s** (crests at
0.065 / 0.135) on a slow residual rise (208.47 → 208.65). Every trailing window
≤ 0.15 still fails G4 (z = 1.15–1.58). σ_inst ≈ 0.1 K, τ_int ≈ 5–10 ms →
projection for the plan's Tier-2 window [0.1, 0.3]: σ_J ≈ 0.028 K →
5√2σ_J ≈ 0.2 K, i.e. an order of magnitude above the ΔJ differences measured
here — identifiability at that window length is **not guaranteed** and must be
re-measured, not assumed.

### Secondary evidence: pre-arrival plateau FD, window [0.01, 0.03]

Recorded only as a consistency check of the FD machinery (the outlet has not
yet seen the developed plume — not a substitute for the planned sensitivity):

- σ_inst = 2.151e-04 K, τ_int = 0.565 ms, σ_J = 5.111e-05 K; half-diff
  6.0e-05 → literal G4 fails but calibrated z = 0.59 (stationary within noise);
  SNR gate passes for both h (threshold 3.61e-04).
- slopes: **+0.028061** (h=0.01) vs **+0.028074** (h=0.02) per α̂ — dev
  **0.048 %**, same sign; curvatures −0.0036 / −0.0098 K (≈ 0);
  dJ/dα ≈ **+2.91e-04 K per m/s**.

### Verdict and next steps

- Tier-1 executed and recorded per plan rules: **G1/G2/G3/G7 pass** (platform
  is bit-reproducible, objective plumbing exact, runs healthy);
  **G4/G5/G6 fail** → per plan §6/§8 the sensitivity in [0.02, 0.06] is
  reported as **not identifiable at this window length** (no gradient claim).
- Escalation path per plan ("move t0 later / widen the window"): Tier-2
  window [0.1, 0.3] at endTime 0.3 (24000 steps, `nStepsFrac` 0.6667 still
  exact). **Deferred by user decision 2026-10-01** (stop at Tier-1).

### Method notes

- Every run logs the per-step objective (`meanTOutlet: <inst> average: <run>`
  lines, 4800 of them) → any trailing window's J is computable offline for
  every run; validated against `evalFunctions` to 5.7e-13. The extra
  `surfaceFieldValue` FO in `d3_j0s` writes a single appending `.dat`
  (needs `name <patch>;` — OF-2506 rejects the older `patch`/`patches`
  selection keys) and matches the log series to text precision (5e-8).
- G4's literal `< 1σ_J` threshold is conservative: for stationary data the
  half-difference has σ = 2σ_J, so ~68 % of stationary windows fail it; the
  analysis script therefore also reports z = halfdiff/(2σ_J).
- Adjoint status (corrects the prompt pack): `libDASolverADR.so` (Sep 21)
  **does** contain 28 `DAContrailFoam` symbols but predates the D2/Prompt-2
  edits; `libDASolverADF.so` (Sep 15) has **none**. Adjoint functionality
  remains unavailable until a joint ADR/ADF rebuild (still gated).

## Phase D4 — ADR/ADF rebuild + gates G1–G5 (in progress, session 2026-10-01)

Execution of `DAContrailFoam_phaseD4_adr_audit.md`. User decisions: Tier-2
[D3] deferred ("costly on my local laptop; go back to Tier-2 later");
option **B** (source patch, below) approved over patchelf-shadow option A
on robustness. **No "ADR support" claim anywhere until gate G5 passes.**

### D4 rebuild (COMPLETE)

- Runner: `work/scripts/build_ad_adf.sh` (→ `/workspace/scripts/build_ad_adf.sh`),
  launched `docker exec -d dafoam-dev bash --noprofile --norc /workspace/scripts/build_ad_adf.sh`.
  OOM-safe: `/tmp/wmake-shim/wmake` rewrites bare `-j` → `-j2` (container has
  10 cores but cgroup `memory.max` = 7.25 GiB; AD `cc1plus` gets OOM-killed at
  10-way parallelism), `pin_env()` re-fronts the shim after every env source,
  `set -e` disabled while sourcing OF env (`pop_var_context`), real status via
  exit codes (upstream `Allmake` `&&`-chain masks failures and its wrapper
  status file lies).
- Result: status `/workspace/logs/dafoam-Allmake-status.txt` =
  `SUCCEEDED time=2026-10-01T11:15:01+00:00`, log
  `dafoam-Allmake-20261001-105945.log`, includes `pip install .`
  (`Successfully installed dafoam-5.1.1`). `dafoam/libs/ADR/` and
  `dafoam/libs/ADF/` both fresh; this supersedes the D3 "stale ADR/ADF"
  note (28-symbol finding) — `libDASolver{ADR,ADF}.so` now postdate D2.

### G1 debugging chain (all fixes recorded, in order)

| # | Failure | Fix |
|---|---------|-----|
| 1 | binary `polyMesh/points` unreadable under codi scalar (`binaryBlock` IO error) | re-`blockMesh` with `writeFormat ascii` |
| 2 | `::sqrt(codi)` etc. hard errors (ADL can't find OF functions for codi types) | unqualified `sqrt/tanh/atan2/sin/cos` in the 21 coded-BC dicts (numerically identical in double mode) |
| 3 | `ld: cannot find -lOpenFOAM -lfiniteVolume -lmeshTools` | root cause + option-B patch (below) |
| 4 | `FOAM FATAL IO: Unescaped '\n' while reading string` (`0/T`) | OF dict strings **cannot contain real newlines**; `\n` escapes are kept **literal** (NOT decoded) → dict-only "second `LIB_LIBS` assignment" trick is impossible; `codeLibs` line removed from all 21 files |
| 5 | bash executed `d4_g1_adr_init.py` as a shell script ("File name too long", `syntax error near '('`) | an env script in the `loadDAFoam`/`AD bashrc` source chain iterates `"$@"` and dot-sources file-like args → `d4_adr_env.sh` now saves `CMD=("$@")`, `set --` before sourcing, `exec "${CMD[@]}"` (also fixed: `exec "$@"` after `set --` was a silent no-op — launch died with exit 0 and empty log) |
| 6 | ADR primal diverges: `DILUPBiCGStab: Solving for h` → `Final residual 1.59e+11` on iteration 1, then `solution singularity` on `h`/`k`/`U` (assembly is fine: `initRes` matches original to 16 digits; `smoothSolver`/`GAMG`/`diagonal` unaffected; identity preconditioner also fails) | **root cause #6, below**: OF's generic `sumProd` uses the `scalarProduct`/`&&` ("dotdot") operator, wrong for scalar-like codi types → `gSumProd` garbage → BiCG/CG `alpha/omega` wrong. Fixed by a scalar `sumProd` specialization in `scalarField.{H,C}` + `libOpenFOAM` rebuild |
| 7 | `wmake libso -j2` → `wmake error: could not change to directory '-j2'` | wmake usage is `wmake [OPTION] [dir]` → **`wmake -j2 libso`** (the `-j`→`-j2` shim is only for the bare `-j` form) |
| 8 | final link fails `/bin/ld: cannot find -lPstream` (`src/OpenFOAM/Make/options:27`; the image ships only `libPstreamADR.so`/`libPstreamADF.so`, in `lib/dummy` + `lib/sys-openmpi`, with no SONAME) | `Make/options` → `-lPstream$(WM_AD_MODE)` (backup `Make/options.orig-d4`); NEEDED comes out `libPstreamADR.so` as before. Note: touching `Make/options` makes wmake regenerate the makefiles → **all 601 objects recompile** (~18 min), not just a relink |
| 9 | union `PYDAFOAM` init path breaks under ADR (double `initSolver()` fatals, pyofm fatals *after* ADR `initSolver`, `setPrimal*` take the double-first path) | official `scripts/d4_g1_adr_init.py` now patches `_initSolver` (construct both bindings, run `_readMeshInfo` **before** ADR `initSolver`, skip double `initSolver`), caches `_readMeshInfo`, and routes both `setPrimal*` to `solverAD` |

### Root cause #6 — OF's generic `sumProd` uses the `&&` (dotdot) operator

The failure was **not** the environment (an ADR-only process that never loads
the original solver still failed identically), not the assembly (after one
`GaussSeidel`-solved step 1 the ADR `p` `initRes` matched the original run to
16 digits), and not DILU specifically (`preconditioner none` failed too).
`GaussSeidel`/`GAMG`/`diagonal` were unaffected because they only need
`gSumMag`/`gSumSqr` (unary) and direct coefficient loops.

The culprit is `src/OpenFOAM/fields/Fields/Field/FieldFunctions.C:488`:

```cpp
TFOR_ALL_S_OP_F_OP_F(resultType, result, +=, Type, f1, &&, Type, f2)
```

which expands to `result += f1[i] && f2[i]` — OF's `scalarProduct` "dotdot"
operator, correct for `Vector`/`Tensor` (where `&&` *is* the dot product) but
a **logical-AND count** for scalar-like types. `double`/`float` mask this with
`*` specializations in `scalarField.{H,C}`; codi `ActiveType` (= `scalar` when
`WM_AD_MODE=ADR/ADF`) had none, so `gSumProd` returned "how many nonzero
pairs", not the inner product → BiCGStab's `alpha = gSumProd(rA,rAyA)/...` ≈ 1
and `omega = gSumProd(tA,sA)/gSumSqr(tA)` garbage → divergence. The same
`&&`-based line exists in v2506, so it is a latent upstream OF bug that only
codi exposes.

**Fix** (option B, backups `scalarField.H.orig-d4` / `scalarField.C.orig-d4`):

- `scalarField.H` (after the `double` declaration): declaration
  `template<> scalar sumProd(const UList<scalar>& f1, const UList<scalar>& f2);`
- `scalarField.C` (after the `double` definition): the definition with `*`
  instead of `&&`, mirroring OF's own `double`/`float` specializations.

A **declaration in the header is mandatory**: `Field.H` includes
`FieldFunctions.C`, so without it every consumer TU inlines the broken generic
(confirmed: pre-fix `libOpenFOAMADR.so` exported no scalar-ActiveType
`sumProd` and no undefined reference either). After the rebuild,
`scalarField.o` defines it (`T Foam::sumProd<codi::ActiveType<...>>(...)`)
and `PBiCGStab.o` carries it as `U` → it calls the specialization.

`scalarProduct<arg1,arg2>::type` is `pTraits<arg1>::cmptType`
(`primitives/VectorSpace/products.H:153`) → for `scalar` it is `scalar`, so
the specialization signature matches the primary template.

Consumers of `gSumProd`/`sumProd` for scalars are only
`{FPCG,PBiCG,PBiCGStab,PCG,PPCG}.C`, `PBiCCCG.C` and `scalarField`/
`complexField` — all inside `libOpenFOAM`, so only that lib had to be
rebuilt (per mode: ADR done, ADF pending).

### Root cause #3 and the approved fix (option B)

OF-AD's `codedFixedValue::prepare()` (and identical code in `codedMixed`)
hardcodes unsuffixed lib names in the generated dynamicCode `Make/options`:

```cpp
// src/finiteVolume/fields/fvPatchFields/derived/codedFixedValue/
//   codedFixedValueFvPatchField.C  (~line 117, same in codedMixed)
      + "\n\nLIB_LIBS = \\\n"
        "    -lOpenFOAM \\\n"          // → patched to -lOpenFOAM$(WM_AD_MODE)
        "    -lfiniteVolume \\\n"       // → -lfiniteVolume$(WM_AD_MODE)
        "    -lmeshTools \\\n"          // → -lmeshTools$(WM_AD_MODE)
      + context.libs()
```

Unsuffixed libs do not exist under `WM_AD_MODE=ADR/ADF` (the AD platform only
has `lib*ADR.so` / `lib*ADF.so`). Patch applied 2026-10-01 to **both** files;
`codeLibs` dict entries (the rejected append-only workaround) removed from the
21 `0/` files of `/tmp/stage3/d4_g1`. `$(WM_AD_MODE)` make-expansion is the
proven pattern (used throughout DAFoam's own `Make/options`). Rejected
alternatives: patchelf SONAME-shadow libs (stale-copy trap, mandatory
LD-order flip = silent ABI hazard), symlink/no-SONAME (mixed process would
bind the original double lib → dlopen undefined symbols).

Known gap (not patched, not needed by this case): `codeStream.C` in
OpenFOAM core has the same hardcode — needed only if a case uses
`#codeStream` (grep count in `d4_g1/0/` = 0).

### G1 PASSED (2026-10-02)

```text
G1: case=/tmp/stage3/d4_g1
G1: ADR initSolver done
G1: INIT_OK elapsed=77.3s
G1: ADR_PRIMAL primalFail=0 elapsed=216.4s
G1: PASS (ADR init + primal + coded-BC compile OK)
EXIT=0
```

Evidence: log `/workspace/logs/g1_gate_official.log`; pristine case
(`/tmp/stage3/d4_g1`, **unmodified DILU/PBiCGStab config** — no solver
downgrades), fresh `dynamicCode` (9 coded BCs + `turbulenceScales_dynamicKEqn`
compiled and linked, `fieldAverage1` ran), 100 steps, **zero**
`singularity`/`nan` solves. Baselines: original mode `primalFail=0` 40 s;
ADR pre-fix diverged on step 1; the `smoothSolver` fallback
(`/tmp/stage3/d4_g1_diag10`, 397 s) is **no longer needed**.

Verified-then-discarded candidates: downgrading all five `PBiCGStab`/`DILU`
entries to `smoothSolver` (passes but deviates from the validated D3 config);
LD_PRELOAD shim (impossible — the pre-fix scalar `sumProd` was fully
inlined, no symbol to override).

### G2 PASSED (2026-10-02) — residual parity, original vs ADR

Cross-process (per the audit: no mixed-process comparison), same case dir and
same states, two envs:

| | original | ADR |
|---|---|---|
| env | `scripts/d4_orig_env.sh` (loadDAFoam only, asserts `WM_AD_MODE` empty) | `D4_AD_MODE=ADR d4_adr_env.sh` |
| script | `scripts/d4_g2_residual_parity.py orig <case>` | same with `adr` |
| log | `logs/g2_residual_orig.log` | `logs/g2_residual_adr.log` |
| result | `n=950516 norm2=5205369.2529319227 absmax=1400533.6289946034` | **identical** |

Verdict (`scripts/d4_g2_compare.py`, log `logs/g2_compare.log`):

```text
G2: max_abs_diff=0
G2: elements differing at all: 0 / 950516
G2: PASS (max elementwise rel 0 <= 1e-10)
EXIT=0
```

Bitwise equality across all 950 516 residual entries — stronger than the
required ≤1e-10. Uses `getResiduals(double*)` (which calls
`updateStateBoundaryConditions()` + `calcResiduals()` itself and returns
full-precision doubles) instead of `calcPrimalResidualStatistics("print")`,
whose `Info` stream only carries 6 significant digits. `dynamicCode/` was
cleared before each run, so the coded BCs compiled fresh against v2506
(`linux64GccDPInt32Opt`, empty `WM_AD_MODE`) and against ADR — confirming the
`$(WM_AD_MODE)` patch works in both modes.

### G3 (2026-10-02) — checkpoint/replay integrity

Case `/tmp/stage3/d4_g3` (copy of `d4_g1`) changed two ways, both required:

- `writeControl timeStep; writeInterval 1` — `writeAdjStates` is gated by
  `runTime.writeTime()` (`DAContrailFoam.C:190-192`), so the inherited
  `adjustableRunTime`/`writeInterval 0.00125` wrote only the final time dir
  and there was no checkpoint set at all.
- `writePrecision 16` (case ships 10) — the audit's ≤1e-12 round-trip bar
  cannot be met at 10 significant digits.

Script `scripts/d4_g3_replay.py`, ADR env, `unsteadyAdjoint.mode =
timeAccurate`, `function.meanTOutlet` with `nStepsFrac 0.5`, log
`logs/g3_replay.log`:

| check | result |
|---|---|
| primal | `primalFail=0`, n=950516 |
| checkpoint set | **101** dirs, `0` … `0.00125` = all 100 steps |
| `getdFScaling` | nonzero at **50..99** (50 of 100), contiguous, all `0.02 = 1/50`; every out-of-window index exactly `0` → 1/N in-window, 0 out-of-window ✓ |
| `primalFinalTimeIndex_` | valid — Prompt-2's fatal did **not** fire ✓ |
| replayed states vary | max step-to-step change > 0 → per-step data really is loaded ✓ |
| final round trip | abs `5.4569682106375694e-12`, scale `16273.161054400569` → **rel `3.3533547614966317e-16`** |

The single `[FAIL]` was my own check comparing the *absolute* round trip to
1e-12; the gate's bar is on the states, and relative it is machine epsilon.
The script now normalizes by `max|primal_final|` and prints `G3: PASS` on a
re-run — the numbers above are the evidence.

**Finding — case-setup gap, must be fixed before G4/G5:** `0/` contains no
`phi`. `readStateVars` `MUST_READ`s every state, and the first three steps
need time levels at/before t=0 that only exist in the hand-provided `0/`
folder (the write sits *inside* the time loop, so the solver never writes
time 0). Observed fatal: `cannot find file ".../0/phi"`. DAFoam's own
`mphys` `readZeroFields` path calls `readStateVars(0.0, deltaT)` before the
primal (`mphys_dafoam.py:1487`), so the unsteady adjoint **cannot run on
this case until every state field (at least `phi`) is in `0/`**. G3
therefore replays the 98 checkpoints whose `t-2·dt` is also a checkpoint
(t ≥ 3·dt); the skip is printed by the script.

### Current state (2026-10-02 ~14:10 UTC)

- **Both AD modes carry the sumProd fix**: `libOpenFOAMADR.so` (24 933 632 B,
  NEEDED `libPstreamADR.so`) and `libOpenFOAMADF.so` (19 944 920 B, NEEDED
  `libPstreamADF.so`), each verified to define exactly one
  `sumProd<codi::ActiveType>` and leave none undefined. Pre-fix copies
  backed up in `/workspace/logs/backups/`
  (`libOpenFOAM{ADR,ADF}.so.pre-sumProd-fix-*`).
- Patched files (backups `*.orig-d4`): `scalarField.{H,C}`,
  `src/OpenFOAM/Make/options`, plus the earlier `codedFixedValue`/
  `codedMixed`, `finiteVolume/Make/options` and the two functionObject
  `$(WM_AD_MODE)` patches.
- Build logs: `libOpenFOAM-ADR-rebuild.log` (first attempt, failed at
  `-lPstream`), `libOpenFOAM-ADR-relink.log` (EXIT=0),
  `libOpenFOAM-ADF-rebuild.log` (EXIT=0).
- Gates: **G1 PASS, G2 PASS, G3 PASS** (above). Next: put `phi` (and any
  other absent state) into `0/` of the G4/G5 case, then **G4**.
- New gate scripts: `d4_orig_env.sh`, `d4_adr_env.sh`,
  `d4_g1_adr_init.py`, `d4_g2_residual_parity.py`, `d4_g2_compare.py`,
  `d4_g3_replay.py`, `d4_rebuild_libOpenFOAM.sh` (all committed to
  `d4_gate_scripts/` in this repo; the originals live in the host
  `work/scripts/`, which is what `/workspace/scripts` maps to).

### Restart commands

```bash
# 1. G1 (re-runnable; clears dynamicCode first so the compile is exercised)
docker exec dafoam-dev bash --noprofile --norc -c 'rm -rf /tmp/stage3/d4_g1/dynamicCode'
docker exec -d dafoam-dev bash --noprofile --norc -c \
  'cd /workspace && D4_AD_MODE=ADR /workspace/scripts/d4_adr_env.sh \
   bash -c "python /workspace/scripts/d4_g1_adr_init.py /tmp/stage3/d4_g1" \
   > /workspace/logs/g1_gate_official.log 2>&1; echo "EXIT=$?" >> /workspace/logs/g1_gate_official.log'
# expect: G1: PASS + EXIT=0 (see G1 PASSED block above)

# 2. ADF rebuild of libOpenFOAM (needs -j2 BEFORE the target)
cp -p /home/dafoamuser/dafoam/OpenFOAM/OpenFOAM-AD/platforms/linux64GccDPInt32OptADF/lib/libOpenFOAMADF.so \
      /workspace/logs/backups/
D4_AD_MODE=ADF /workspace/scripts/d4_adr_env.sh bash /workspace/scripts/d4_rebuild_libOpenFOAM.sh \
  > /workspace/logs/libOpenFOAM-ADF-rebuild.log 2>&1
# wmake emits the unsuffixed name (Make/files has no AD suffix) -> rename it
# in the same dir:
#   mv .../platforms/linux64GccDPInt32OptADF/lib/libOpenFOAM.so \
#      .../platforms/linux64GccDPInt32OptADF/lib/libOpenFOAMADF.so
# then verify: nm -C --defined-only ... | grep -c "sumProd<codi::ActiveType" == 1

# 3. G3 (re-runnable; primal ~5 min, replay ~15 min — ascii 16-digit
#    checkpoint reads are I/O bound, ~4.8 GB total). Expect the [OK] lines
#    for checkpoints / dFScaling window / in-window 1/N / varying states,
#    round-trip ~3.35e-16 rel, then "G3: PASS" and EXIT=0.
docker exec -d dafoam-dev bash --noprofile --norc -c \
  'rm -rf /tmp/stage3/d4_g3/dynamicCode; cd /workspace && \
   D4_AD_MODE=ADR /workspace/scripts/d4_adr_env.sh \
   bash -c "python /workspace/scripts/d4_g3_replay.py /tmp/stage3/d4_g3" \
   > /workspace/logs/g3_replay.log 2>&1; echo "EXIT=$?" >> /workspace/logs/g3_replay.log'

# watch it (host is zsh — quote everything):
docker exec dafoam-dev bash --noprofile --norc -c \
  'grep -n "G3:\|EXIT=" /workspace/logs/g3_replay.log | tail -20'
```

**Before G4/G5**, the case's `0/` must contain every state field including
`phi` (G3 finding above): `cp <ckpt>/0.00125/phi <case>/0/`, then re-run
G3 and require `0/ missing state fields: []`. `dynamicCode` must be removed
whenever `WM_AD_MODE` changes.

Case `d4_g1` = copy of `d3_base` (both container-local `/tmp/stage3/`) with
`endTime`/`writeInterval` 0.00125 (100 steps), `writeControl
adjustableRunTime`, ascii mesh, unqualified math in the BCs, no `codeLibs`.
If container `/tmp` is lost: `cp -r /tmp/stage3/d3_base /tmp/stage3/d4_g1`
recreates the base, then re-apply those four edits (and `blockMesh` with
ascii `writeFormat`). `/tmp/wmake-shim/wmake` (417-byte `-j`→`-j2` rewrite
wrapper, falls back through `$PATH` for the real wmake) must be recreated
if `/tmp` is wiped.

### Gates remaining

- **G1 PASS (2026-10-02)** — ADR build + coded-BC compile + 100-step ADR
  primal (evidence above).
- **G2 PASS (2026-10-02)** — residual parity, bitwise identical across
  modes: 950 516 residuals, `max_abs_diff = 0` (evidence above).
- **G3 PASS (2026-10-02)** — 101 checkpoints, `getdFScaling` 1/N in-window
  (50..99) and 0 outside, `primalFinalTimeIndex_` valid, final-state round
  trip **3.35e-16 relative** (evidence above). Blocked the first 3 steps'
  oldTimes: `0/` has no `phi` — **case-setup fix required before G4/G5**.
- **G4** derivative plumbing without chaos — (a) short-window `meanTOutlet`
  derivative for an existing differentiable input on a fixed-dt window;
  (b) instrument the §7 non-smooth inventory for activity in that window.
- **G5** FD acceptance — fixed-dt window `endTime 0.012`, objective
  `meanTOutlet`, control = jet velocity scale α, central difference with
  h ∈ {0.005, 0.01, 0.02}, `|dJ/dα|_AD − |dJ/dα|_FD| / |dJ/dα|_FD ≤ 1e-4`
  with D3's R1–R4 gates (audit §9).
- **G6** record — G1–G5 evidence in the runbook, prompt-pack ADR/ADF bullets
  flipped to verified.

Record every gate's evidence here; ADR-support claim only after G5.

### Environment gotchas (verified this session)

- `d4_adr_env.sh` = union env: `loadDAFoam.sh` + `WM_AD_MODE` sed +
  `OpenFOAM-AD/etc/bashrc` + original v2506 lib dirs appended (AD dirs stay
  first; SONAMEs disjoint so both instances' `.so`s resolve in one process).
  `D4_AD_MODE=ADR|ADF` now parameterizes it.
- OF dict strings: real newlines rejected, `\n` NOT decoded → never rely on
  newline embedding in `.dict` values.
- wmake objects live under `OpenFOAM-AD/build/<WM_OPTIONS>/src/...`; keep
  `-j` ≤ 2 (OOM).
- Host is zsh: bare `===` in a command → `== not found`; `pgrep -f` matches
  its own command line (filter with `pgrep -x` / `grep -v grep`).

## Notes

The implementation should stay intentionally narrow and reviewable. The first milestone is not a full DAFoam contrail solver; it is the reduced gas-phase path that matches the validated CASSANDRA no-PBE physics.
