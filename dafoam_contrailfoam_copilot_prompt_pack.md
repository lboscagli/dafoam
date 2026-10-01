# DAContrailFoam Local-LLM Prompt Pack

## Purpose

Use these prompts in Continue, Copilot Chat, or another code assistant while developing the CASSANDRA–DAFoam integration. They are intentionally narrow: request **one reviewable change at a time**, inspect the diff, compile, test, then commit.

This document is a prompt checklist, not permission to accept large automatic edits. The developer remains responsible for reviewing all OpenFOAM/DAFoam changes.

## Current validated baseline

### Repositories

| Repository | Role |
|---|---|
| `cassandra` | Standalone `contrailFoam`, `contrailCellFoam`, `aerosolMicrophysics`, `momentTransport`, cases, scripts, and CASSANDRA documentation. |
| `dafoam` | DAFoam-specific implementation: `DAContrailFoam`, state/residual support, functions, and DAFoam tests. |

### Local Docker paths

```text
CASSANDRA: /workspace/src/cassandra
DAFoam:    /workspace/src/dafoam
```

Load CASSANDRA/OpenFOAM in a clean container shell with:

```bash
source /workspace/src/cassandra/dafoam_env_doc/load_cassandra_dafoam_env.sh
```

### Version and branch baseline

```text
DAFoam version:       v5.1.1
DAFoam baseline:      a941594dc3d0a85b885074f8aabcf2806ad4604d
DAFoam feature branch: feature/contrailfoam-primal-wrapper
OpenFOAM:             v2506
```

### Operational status in this workspace

Progress log (last updated after Phase C validation and the outlet-BC investigation):

**Done and verified:**

- `DAStateInfoContrailFoam` is present and registers the reduced gas-phase state set (`U, p, T, O2, CO2, H2O, nut` per cell, `phi` per face; `N2` reconstructed, not a state).
- `DADynamicKEqn` header/source skeletons are present and registered in the DAFoam runtime-selection build list; **its kResidual dimensions are fixed (Phase A)** so `calcPrimalResidualStatistics` passes (`Total Residual Norm2 12447191.93`).
- `DATurbulenceModel` no longer assumes a hard-coded `RAS` dictionary and now accepts an `LES` turbulenceProperties block when present.
- `DAResidualContrailFoam` implementation exists, is registered, and compiles.
- `DAContrailFoam` wrapper compiles and links in original mode: Stage 1 fixed the two root compile errors (`tmp<convectionScheme>::ptr()` deletion in `createRefsContrail.H`, removed `GeometricField` ctor in `YEqnContrail.H`); Stage 2 rewrote `createFieldsContrail.H` (thermo/species/inert validation, rho/U/phi, pressureControl, turbulence, dpdt/K, MRF), fixed `initSolver()` ordering (pimple control -> fields -> LES-or-RAS turbulence name -> `DATurbulenceModel::New` -> `createAdjoint.H`), and fixed `createRefs` ownership (`thermo.p()`, no null `T`).
- Registration is wired: `src/adjoint/Make/files` entry and the `pyDAFoam.py` `Compressible` registry entry.
- Original-mode build validation passed: log `logs/dafoam-Allmake-20260929-152848.log`, 0 errors (authoritative check; `dafoam-Allmake-status.txt` can be stale).
- PYDAFOAM initialization smoke test passed on the LES reference case (case `0/` generated via `generate_inlet_BCs.py`, mesh built via `blockMesh`): full construction (`initSolver` -> `createFields` -> turbulence -> `createAdjoint`, every registered-object lookup succeeded), exact adjoint state count `9*nCells + nFaces`, and `dpdt`/`K`/`rho` readable. Smoke script: `scripts/smoke_test_dacontrail.py` (original mode only; ADR/ADF bypassed).
- **Stage 3 / Phase B:** `solvePrimal` time loop + CASSANDRA equation parity committed as `6f7a7ca` "Add Stage 3 primal time loop and contrail equation parity" (9 files, +288/−104, incl. `rhoEqnContrail.H`). B5's "preserve the transient CASSANDRA gas-phase PIMPLE sequence exactly" is satisfied; `fvSolution` has `consistent no`, so the `pcEqn` branch is inert in both solvers.
- **Phase C executed (engineering-level pass):** step gates (step 1 bit-exact except FP-level k source; step 3 `U` 5.75e-6 / `k` 1.35e-2 / `nut` 1.88 — chaotic Germano amplification); matrix bisect proved the k seed is FP-level source-vector differences only (kEqn diag and all other dumped inputs bit-identical); twin runs `da_long`/`of_long` to t=0.076 with a common write at t=0.06: inlet patch values bit-identical, means agree to ≤3.5% (`UMean`) / ≤1.0% (`TMean`, `rhoMean`), figures exported to `work/stage3_figures/` (47 vs 50 PNGs, visually near-identical per user).
- **Outlet-BC investigation closed with no defect:** `0/` byte-identical; divergence grows monotonically inlet (FP-level) → outlet; worst diffs in plume interior; outlet biases physically consistent with plume-arrival differences (see runbook "Phase C validation findings"). VTK suffix difference (`da_long_1` vs `of_long_1685`) is write-metadata only: OF's full write stores `<time>/uniform/time` `index 1685`, DA's reduceIO write omits it so `foamToVTK` falls back to directory enumeration; both dirs are t=0.06 (`.vtm.series`).
- The OpenFOAM-v2512 tree was used only as a reference for the dynamicKEqn formulas; it was not compiled and is not part of the target build.
- **Phase D2 complete (commit `f33eab6`):** `functionTimeSteps_` bounds guard in `DASolver::calcAllFunctions` (grow-on-demand + `listIndex >= 0` skip); original-mode rebuild clean (`logs/dafoam-Allmake-20261001-062527.log`, `*** Build Successful! ***`); `stage3_primal_compare.py --func` injects `meanTOutlet` and prints `evalFunctions`; 0.012 s validation vs `fieldAverage1` → `evalFunctions` 208.0155615012 K vs `TMean` outlet areaAverage 208.0159243 K, **difference −3.628e-4 K = −1.744e-6 rel, entirely dt-weighting** (`DATimeOpAverage` unweighted, fieldAverage dt-weighted; both sample `t1..tN`, no t=0 sample); guard-trigger run (60-slot list, 120 steps) survived with exact readback. Full record: runbook "Phase D2".

**Remaining / not yet complete:**

- Phase C formal closure (optional): 0.3 s runs were stopped at t≈0.076; re-run for a full-length final/time-mean comparison if a formal sign-off beyond the current evidence is wanted.
- Phase D: **D1/D2/D3/D4 all delivered.** D1 design
  (`DAContrailFoam_phaseD1_mean_objective_design.md`); **D2 complete (commit
  `f33eab6`)** — guard + `--func` + validation, dt-weighting difference
  quantified (−1.74e-6 rel); **Prompt 2 complete (commit `fef5cc5`)** —
  fail-fast errors for primal-not-run / unknown function / empty window /
  zero patch area, verified by five failure-case runs (runbook "Phase D2
  hardening"); **D3 plan**
  (`DAContrailFoam_phaseD3_sensitivity_plan.md`): inlet jet velocity scale α,
  fixed-dt FD protocol with repeatability gates R1–R4 and acceptance G1–G7,
  costed two-tier run matrix; **D4 audit**
  (`DAContrailFoam_phaseD4_adr_audit.md`): non-smooth inventory, coded-BC
  blocker, checkpoint/replay + time-average requirements, gates G1–G6.
  Remaining work is *execution*, not design: run the D3 matrix; on explicit
  request rebuild ADR/ADF starting at D4 gate G1 with acceptance at G5.
  Session state + resume prompts: `DAContrailFoam_session_handover.md`.
- ADR/ADF: AD libs are stale (`libDASolverADR.so` Sep 21, `libDASolverADF.so` Sep 15, no `DAContrailFoam`); the D4 audit is written (`DAContrailFoam_phaseD4_adr_audit.md`) — rebuild only when explicitly requested, then gates G1→G5 apply (no "ADR support" claim before G5).
- `reduceIO` intermediate writes omit `N2`, `nut`, `alphat`, `rho` (and `uniform/time`); comparison scripts must not expect them in intermediate time dirs.
- `DAField.C:1137` hard-codes `turbDict.subDict("RAS")` — dormant because default `primalBC = {}` (fires only with `useWallFunction`), but must be generalized before wall-function LES use.
- aerosol/PBE physics (Phase E, intentionally deferred)

**Working tree:** this changeset commits the Phase A–D2 record (`DAContrailFoam_runbook.md`, this progress update, `DAContrailFoam_phaseD1_mean_objective_design.md`, `DAContrailFoam_session_handover.md`) together with the D2 bounds guard in `src/adjoint/DASolver/DASolver.C`. The Stage 3 code itself was committed earlier as `6f7a7ca`. Nothing further is uncommitted; further commits only on explicit request.

### Completed DAFoam changes

The following commits are compiled in original mode (local branch; push only on explicit request):

```text
759e6d9 Allow multicomponent thermophysical dictionaries
8046cff Add contrail gas-phase state information
6f7a7ca Add Stage 3 primal time loop and contrail equation parity
f33eab6 Add mean objective bounds guard and phase D2 record
```

`DAStateInfoContrailFoam` currently registers the independent reduced-model gas states:

```text
U, p, T, phi, O2, CO2, H2O, and runtime-corrected turbulence model states
```

`N2` is not an independent state: it is reconstructed as the inert species.

> **LES status:** the **primal** `dynamicKEqn` LES path is now live: `DADynamicKEqn` is registered, its kResidual dimensions were fixed in Phase A, and `DAContrailFoam` has run the dynamicKEqn LES reference primal (Phase C twin runs). This is primal evidence only — it is **not** evidence that reverse/forward AD through `dynamicKEqn` works; that remains gated on the Phase D4 ADR enablement audit.

### First physical target

The intended scientific target is the **gas-phase-only LES** configuration:

```text
aerosolActive      false;
massCouplingActive false;
heatCouplingActive false;
nAerosolModes      0;
simulationType     LES;
LESModel           dynamicKEqn;
```

LES reference case:

```text
/workspace/src/cassandra/run/contrailFoam_nasa_pal_HX1/
nasa_pal_figure14_randomU
```

or a reduced, copied development case derived from it.

The target is transient, compressible, multicomponent, and uses:

```text
psiReactionThermo
multiComponentMixture
species: O2, N2, CO2, H2O
inert species: N2
```

### Mandatory LES prerequisite

The DAFoam fork currently provides turbulence-adjoint classes for RANS models such as `kEpsilon`, but not a validated `dynamicKEqn` LES implementation. The implementation order for the LES target is:

```text
[Completed in local workspace]
1. Audit the OpenFOAM dynamicKEqn structure and DAFoam turbulence interface.
2. Add the DAFoam dynamicKEqn model scaffold and register it in the build.
3. Generalize DATurbulenceModel to accept LES dictionaries without assuming RAS.
4. Validate the adjoint-model code by compiling original mode.        <-- done
   (build log dafoam-Allmake-20260929-152848.log: 0 errors, DAContrailFoam
   wrapper included; see Operational status for the Stage 1+2 patch)
5. Add and validate the multicomponent gas residual path.              <-- partially
   (DAResidualContrailFoam implemented and compiling; evaluation still
   aborts on a dimensions mismatch, see Operational status)

[Still required]
6. Run a minimal LES primal smoke test using the CASSANDRA gas-phase
   reference case (initialization smoke test passed; the primal itself
   still needs the Stage 3 time-loop fix before it can be run).
7. Add unsteady, time-averaged LES observations.
8. Only then plan ADR/adjoint support.
```

Do not falsely use the `dummy` model, a RAS dictionary, or a `kEpsilon` residual as a substitute for `dynamicKEqn` while claiming LES support.

## Global instructions to prepend to every prompt

Paste this instruction before each task-specific prompt:

```text
Work only in the DAFoam fork unless the request explicitly says otherwise.
Do not edit generated lnInclude/, Make/linux*/, build/, platforms/, case-output,
or processor directories. Do not edit CASSANDRA source unless explicitly requested.
Do not change more files than necessary. First explain the proposed change and list
all affected files; do not apply edits until I explicitly request a patch.

The current implementation target is original-mode, PBE-disabled, multicomponent
LES with dynamicKEqn. Do not add aerosol/PBE, ADR, ADF, objectives, or
data-assimilation code unless explicitly requested. Do not claim or emulate
dynamicKEqn support using a RAS model, a dummy turbulence model, or a RAS
dictionary. First establish the dynamicKEqn turbulence-model prerequisite.

The CASSANDRA species equation solves O2, CO2, and H2O independently, with N2
reconstructed as the inert species. Do not represent the PtrList of species as a
single state named Y. Do not solve equations inside DAResidual::calcResiduals;
residual code must evaluate residual fields for the supplied state.
```

## Mandatory review and validation loop

After every proposed patch:

```bash
git diff
git diff --check
```

After every successful reviewable build:

```bash
git status
git add <only intended files>
git diff --cached
git commit -m "<imperative message>"
git push
```

For initial implementation increments, compile original mode only:

```bash
docker exec -d \
  -e COMPILE_DAFOAM_NOAD=1 \
  --workdir /workspace/src/dafoam \
  dafoam-dev \
  bash --noprofile --norc /workspace/scripts/build_dafoam_source.sh
```

Do not enable ADR/ADF until CASSANDRA custom code and the complete residual path are deliberately ported and derivative-checked.

---

# Phase 0 — dynamicKEqn LES turbulence prerequisite

This phase is mandatory for the LES target. DAFoam v5.1.1 does not currently contain a fully validated `DATurbulenceModel` child for `dynamicKEqn`.

## Status snapshot

Completed in the current workspace:

- the DAFoam `dynamicKEqn` model declaration/registration was created
- the `DATurbulenceModel` base was generalized to read either `RAS` or `LES` subdicts
- the OpenFOAM-v2512 dynamicKEqn formulas were used as implementation guidance
- original-mode build validation passed (including the `DAContrailFoam` wrapper,
  Stage 1+2 patch; build log `dafoam-Allmake-20260929-152848.log`, 0 errors)
- PYDAFOAM initialization smoke test passed on the LES reference case
  (`scripts/smoke_test_dacontrail.py`, original mode)
- **running the LES primal itself** (Stage 3 time loop, commit `6f7a7ca`)
- **residual-evaluation validation** (Phase A kRes dimensions fix;
  `--residuals` passes)
- **standalone-versus-DAFoam LES primal comparison** (Phase C: gates + twin
  0.06 s runs + figures; see Operational status)

Not yet complete:

- any ADR/adjoint work through the LES model layer (Phase D4 audit gate)
- Phase C formal full-length (0.3 s) closure, if wanted beyond current evidence

## Prompt L1 — audit the OpenFOAM dynamicKEqn model and DAFoam turbulence interface

```text
Do not edit files. Compare OpenFOAM v2506 dynamicKEqn source and headers with
DAkEpsilon, DATurbulenceModel, DAModel, DAStateInfoRhoPimpleFoam, and the
CASSANDRA LES reference case.

Identify the actual dynamicKEqn state fields, residual equations, required
thermodynamic/turbulence fields, runtime dictionary entries, and all methods
that a DADynamicKEqn class must implement. State whether the model is
compressible-compatible in this OpenFOAM installation. Do not guess APIs and
do not write code.
```

**Review gate:** List exact OpenFOAM source/header paths and all affected DAFoam classes before proposing a patch.

## Prompt L2 — LES model-state and residual design

```text
Do not edit files. Propose the smallest complete DAFoam design for
DADynamicKEqn support. It must state:
- the independent LES model state fields;
- their residual field names;
- their connectivity to U, p, T, phi, O2, CO2, H2O, and thermodynamic fields;
- how turbulent viscosity and turbulent thermal diffusivity are updated;
- whether existing DATurbulenceModel RAS assumptions must be generalized.

Do not reuse a kEpsilon or dummy model as an approximation. Do not write code.
```

**Review gate:** The proposal must preserve the OpenFOAM `dynamicKEqn` primal model, not merely compile.

## Prompt L3 — minimal dynamicKEqn class declaration and registration

```text
Create only the DADynamicKEqn header and runtime-selection/build registration
required by the approved design. Do not change DAContrailFoam yet. Do not
implement residual equations, species transport, ADR, ADF, PBE, or objectives.

Before editing, list every file and show the intended patch. Do not touch
lnInclude or Make/linux directories.
```

**Build gate:** Original mode only.

## Prompt L4 — dynamicKEqn original-mode implementation and primal validation

```text
Implement only the approved original-mode DADynamicKEqn methods needed for
DAContrailFoam primal execution. Reuse exact OpenFOAM v2506 dynamicKEqn APIs.
Do not claim adjoint support, do not enable ADR/ADF, and do not add PBE.

Provide a standalone-versus-DAFoam LES-primal validation plan that compares
the turbulence fields, nut, alphat, U, p, and T on a short copied case.
```

**Build gate:** Original mode only. Run an isolated LES primal smoke test.

---

# Phase A — Reduced gas-phase residual implementation

**Status:** **Passed.** Prompts A1–A7 executed in earlier sessions; the Phase A blocker (kResidual dimensions `[1 2 -3]` vs `[1 -1 -3]`) was fixed in `DADynamicKEqn.C` (`rhoDimensions()*k_.dimensions()/dimTime`); smoke test and `calcPrimalResidualStatistics`/`--residuals` pass (`Total Residual Norm2 12447191.93`).

Begin this phase only after the `dynamicKEqn` prerequisite supports the LES primal path.

## Prompt A1 — inspect equivalent DAFoam residual patterns

```text
Do not edit files. Compare DAResidualRhoPimpleFoam, DAScalarTransportFoam,
DAStateInfoContrailFoam, the completed DADynamicKEqn model support, and
CASSANDRA's YEqn_aero.H, hEqn_aero.H, UEqn.H, and pEqn.H for the
PBE-disabled LES case.

Propose the smallest complete design for DAResidualContrailFoam. It must include
URes, pRes, TRes, phiRes, O2Res, CO2Res, and H2ORes. Explain exactly which
OpenFOAM fields are derived rather than independent states, how N2 is
reconstructed, and which existing DAFoam code can be reused unchanged.
Do not write code.
```

**Review gate:** Confirm that the proposal does not use the constant-property `DAResidualRhoPimpleFoam::updateIntermediateVariables()` formulation.

## Prompt A2 — residual class header only

```text
Create only DAResidualContrailFoam.H. Do not modify Make/files, pyDAFoam.py,
or create a solver wrapper yet.

The class must declare references/residual fields for U, p, T, phi, O2, CO2,
and H2O; references to psiReactionThermo-compatible thermodynamics, rho, psi,
he, alphat, dpdt, K, compressible turbulence, DADynamicKEqn support, and the
required virtual methods. Use DAResidualRhoPimpleFoam.H as a structural
reference, but do not inherit its constant-property assumptions.

State explicitly which header includes are required for psiReactionThermo,
basicSpecieMixture, multivariate convection, MRF, fvOptions, and pressure
correction. Do not apply the patch until I approve it.
```

**Review gate:** Ensure `O2Res`, `CO2Res`, and `H2ORes` will have names expected by `DAField`.

## Prompt A3 — residual constructor and clear method

```text
Implement only the DAResidualContrailFoam constructor and clear() method.
Register residual fields with setResidualClassMember macros for U, p, T, phi,
O2, CO2, and H2O. Use dimensions consistent with the CASSANDRA finite-volume
species residual. Obtain psiReactionThermo from the mesh registry and obtain
composition/Y species fields through its actual API.

Do not implement calcResiduals, updateIntermediateVariables, boundary
correction, or preconditioner assembly yet. Do not add code that solves an
equation. Explain the field lookup choices before editing.
```

**Build gate:** Add only the residual source to `src/adjoint/Make/files`, build original mode, and commit if successful.

## Prompt A4 — multicomponent intermediate-variable update

```text
Do not edit files yet. Inspect OpenFOAM v2506 psiReactionThermo and
multiComponentMixture APIs used by CASSANDRA. Propose a correct
DAResidualContrailFoam::updateIntermediateVariables() sequence after DAFoam
writes U, p, T, O2, CO2, and H2O from the state vector.

The sequence must reconstruct N2 = 1 - O2 - CO2 - H2O, enforce the same
non-negativity treatment as the standalone solver, update thermodynamic fields
through the actual OpenFOAM thermo model, and update rho, he, psi, K, dpdt,
and relevant boundary conditions. Do not use constant Cp, molecular weight, or
hand-coded Sutherland formulas. Do not write code until I approve it.
```

**Review gate:** The proposal must match OpenFOAM APIs available in v2506, not APIs guessed from other versions.

## Prompt A5 — implement thermo update and boundary correction

```text
Implement the approved multicomponent updateIntermediateVariables() and
correctBoundaryConditions() in DAResidualContrailFoam only. Boundary correction
must include U, p, T, O2, CO2, H2O, and reconstructed N2 as appropriate.
Do not add aerosol/PBE logic. Do not call YiEqn.solve() or modify the solution
inside calcResiduals.

Before applying, show the exact patch and explain how the implementation avoids
the legacy DAResidual::updateThermoVars() path.
```

**Build gate:** Original-mode build. Then run a small source-level inspection or unit initialization test if available.

## Prompt A6 — species residual evaluation

```text
Implement only the O2, CO2, and H2O residual evaluation in
DAResidualContrailFoam::calcResiduals(). Mirror CASSANDRA's YEqn_aero.H with
PBE disabled: ddt(rho,Yi), multivariate convection using phi, and
laplacian(muEff/Sc,Yi). Do not include vaporSink or fvOptions source terms
unless they are demonstrably required in the reduced reference case.

The method must evaluate O2Res, CO2Res, and H2ORes only. It must not call
solve(), clamp any state, reconstruct N2, or advance time. Show the equation
mapping and the patch before editing.
```

**Review gate:** Verify the diffusion coefficient is dynamic viscosity `muEff/Sc`, not kinematic viscosity `nuEff/Sc`.

## Prompt A7 — flow, energy, pressure, and flux residual evaluation

```text
Implement the remaining URes, TRes, pRes, and phiRes portions of
DAResidualContrailFoam::calcResiduals(), matching the PBE-disabled CASSANDRA
UEqn.H, hEqn_aero.H, and pEqn.H rather than blindly copying
DAResidualRhoPimpleFoam. Include the dynamicKEqn-compatible turbulence stress
and thermal-diffusion paths, actual MRF/fvOptions behavior required by the
reference case, and the correct pressure equation branch.

Do not implement an approximate preconditioner beyond the inherited safe
behavior unless it is separately designed. Do not add PBE, aerosol, or
microphysics source terms.
```

**Build gate:** Original mode only. Inspect residual fields at a standalone converged LES state before proceeding.

---

# Phase B — Primal DAContrailFoam wrapper

**Status:** **Complete.** B1–B4 done (file split, equation includes ported from CASSANDRA, field creation per B3 in `createFieldsContrail.H`); B5 done for registration/initialization (Stage 1+2) and for the transient sequence: the Stage 3 rework of `solvePrimal()` (outer time loop, `UEqn`->`YEqn`->`hEqn` order, `rhoEqnContrail`, storePrevIter/CourantNo/setDeltaT/turbulence-correct/write, ddtCorr fix in `pEqnContrail.H`) is committed as `6f7a7ca`. Parity notes: `fvSolution` `consistent no` makes the `pcEqn` branch inert in both solvers; `pEqnContrail.H` ports the non-transonic branch by design (matches `DAResidualContrailFoam`).

## Prompt B1 — design the wrapper file split

```text
Do not edit files. Propose the smallest file split for a gas-phase-only
DAContrailFoam wrapper. It must follow the CASSANDRA no-PBE time loop and use
the completed dynamicKEqn LES support. Identify what should be copied/adapted
into DAFoam and what can be shared safely. Account for the DADynamicKEqn
DAFoam turbulence model. Do not write code.
```

## Prompt B2 — CASSANDRA include/link strategy

```text
Do not edit files. Determine whether the no-PBE DAContrailFoam wrapper requires
CASSANDRA headers or links to libaerosolMicrophysics and libmomentTransport.
Inspect the actual gas-only source includes and Make/options. Propose the
smallest original-mode include/link changes, if any. Do not link unused aerosol
libraries merely because standalone contrailFoam does.
```

## Prompt B3 — field creation implementation

```text
Implement only the DAContrailFoam field-creation include/header needed for the
reduced gas-phase LES case. It must create psiReactionThermo,
basicSpecieMixture, species list, inert species validation, rho, U, p, phi,
pressure control, compressible turbulence, multivariate field table, dpdt, K,
MRF, and fvOptions. It must not create aerosol moments, vaporSink,
latentSource, or microphysics objects.

Use CASSANDRA createFields.H as the physics source of truth. Before editing,
show the exact file list and patch.
```

## Prompt B4 — primal gas-phase equation includes

```text
Implement only the PBE-disabled gas-phase equation includes for DAContrailFoam:
rho equation, momentum equation, active-species transport and inert-species
reconstruction, enthalpy equation, pressure equation, and dynamicKEqn
LES-compatible turbulence correction.

Match CASSANDRA UEqn.H, YEqn_aero.H, hEqn_aero.H, and pEqn.H in the reduced
configuration. Preserve multicomponent behavior and use the same species
names/inert reconstruction. Do not include aerosol transport or PBE calls. Do
not use the DARhoPimpleFoam equation headers unchanged.
```

## Prompt B5 — solver class and registration

```text
Implement the DAContrailFoam class only after the approved field/equation
includes and DAResidualContrailFoam exist. Add runtime selection, Make/files
registration, and the pyDAFoam Compressible registry entry.

In solvePrimal(), preserve the transient CASSANDRA gas-phase PIMPLE sequence
exactly. Keep unsteady input updates, function evaluation, validation, and
reduceIO behavior only where they are compatible with the CASSANDRA order.
Do not enable ADR/ADF and do not add aerosol/PBE code.
```

**Build gate:** Original mode only, then run PYDAFOAM initialization against a copied LES case.

---

# Phase C — Standalone versus DAFoam primal validation

**Status:** **Executed — engineering-level pass (formally partial).** C1/C2 in use: `scripts/stage3_primal_compare.py` (modes `da`/`of`/`cmp`, `k/TKE` + `k/TKEsum` rows) drives both solvers and dumps npz. Evidence: step-1 gates bit-exact except the FP-level k source seed (matrix bisect: kEqn diag + all other dumped inputs bit-identical; only explicit source differs, max 1.56e-8); step 3 `U` 5.75e-6 / `k` 1.35e-2 / `nut` 1.88 (chaotic Germano amplification — use engineering tolerances, not rtol 1e-6); twin runs `da_long`/`of_long` to t≈0.076 with common write at t=0.06: inlet patch values bit-exact, means `UMean` 3.5% / `TMean` 0.8% / `rhoMean` 1.0%, second moments 5–10%, figures in `work/stage3_figures/` (user-verified near-identical). Outlet-BC suspicion investigated and closed with no defect; VTK suffix difference is write-metadata only (see runbook findings). **Optional remaining:** full 0.3 s run + final/time-mean comparison for a formal sign-off. C3 (regression case) not started.

## Prompt C1 — minimal PYDAFOAM driver

```text
Do not edit DAFoam C++ code. Create a minimal Python driver/test for the
CASSANDRA gas-phase LES reference case using solverName DAContrailFoam.
The driver must run primal only, use the existing case directory or a copied
case, and write results to a separate output location. Do not define an
objective or adjoint solve yet.
```

## Prompt C2 — comparison script

```text
Create a comparison script that compares the completed standalone
contrailFoam gas-phase LES reference with the DAContrailFoam primal result at
the final physical time. Compare at least U, p, T, rho, O2, CO2, H2O, and the
relevant dynamicKEqn turbulence fields. Also compare a selected short-window
time average where fieldAverage data are available. Define tolerances and
report maximum and relative L2 differences. Do not modify solver code.
```

## Prompt C3 — DAFoam regression test

```text
Do not edit code yet. Propose a small DAFoam regression test case and test
driver. It must be much smaller than the scientific LES case, remain
multicomponent and PBE-disabled, use dynamicKEqn, and check a small set of
deterministic reference quantities. Explain how sampling uncertainty or
short-window averaging will be handled. Do not add it until the primal wrapper
matches the standalone reference.
```

---

# Phase D — Time-averaged LES objective, then adjoint

**Status:** **D1 done, D2 recon done, implementation pending.** The Phase C gate is satisfied at engineering level (dynamicKEqn LES primal runs and matches standalone within declared tolerances; outlet-BC investigation closed clean). D1 design delivered (`DAContrailFoam_phaseD1_mean_objective_design.md`); D2 recon found the objective machinery already exists — `DAFunctionPatchMean` (dict: `type/source/patches/scale/varName/varType/index`), `DATimeOpAverage`, `getTimeOpRange` trailing `nStepsFrac` window (emulate [t0,t1] with `endTime=t1`, `nStepsFrac=(t1-t0)/endTime`), python `evalFunctions` → `getTimeOpFuncVal` — so D2 = one bounds guard (`DASolver.C:369`, unguarded store into a list sized `round(endTime/deltaT)` at :584) + `function` dict in the driver + short validation run vs the case's `fieldAverage1` (`controlDict` line 377; dt-weighted, unlike `DATimeOpAverage`). ADR/ADF rebuild remains gated on Prompt D4 + explicit request. Work order: D2 (implement + validate) → D3 (sensitivity plan) → D4 (ADR audit) → rebuild + FD validation. Resume from `DAContrailFoam_session_handover.md` §3.

Do not begin this phase until the `dynamicKEqn` LES primal path and the multicomponent gas-phase primal comparison both pass.

## Prompt D1 — design only: mean-observation objective for unsteady LES

```text
Do not write code. Design a DAFoam objective for time-averaged LES
observables. The governing equations and adjoint target remain instantaneous
and unsteady; only the observation operator is averaged over a post-transient
sampling window.

The first objective should support one mean temperature or velocity observable
on a patch, cell zone, sampling plane, or probe-like region. Define discrete
accumulation, reset, normalization, storage, and comparison with measured mean
data. Explain why the objective must compare the mean prediction to the mean
measurement, rather than average instantaneous squared errors.
```

## Prompt D2 — time-averaged LES objective implementation

```text
Implement the smallest objective-function extension required for one
post-transient time-averaged scalar LES observable. Integrate it with
DAFoam's unsteady function evaluation without adding PBE physics. Show all
modified files and explain reset, accumulation, sampling window, finalization,
and restart behavior.
```

**Validation gate:** Compare the computed mean with OpenFOAM fieldAverage or an independently calculated mean from the standalone LES reference.

## Prompt D3 — LES sensitivity and identifiability plan

```text
Do not enable ADR yet. Propose a one-control, one-objective plan for the
gas-phase LES model. Use a scalar control such as ambient H2O mass fraction or
inlet velocity scale and a time-averaged temperature or velocity objective.
Specify a post-transient averaging window, finite-difference perturbation
sizes, repeatability checks, expected turbulent sampling uncertainty, and
acceptance criteria. Do not write code until the plan is approved.
```

## Prompt D4 — ADR enablement audit for dynamicKEqn LES

```text
Do not edit code. Audit all DAContrailFoam, DADynamicKEqn,
DAResidualContrailFoam, multicomponent thermo, CASSANDRA equation headers,
and external libraries for ADR compatibility. Identify every required
original-to-ADR build change and every non-smooth operation. Include unsteady
checkpoint/replay and time-averaged objective requirements.

Produce a checklist and build/test plan. Do not claim ADR support before
finite-difference validation of one short-window mean objective.
```

---

# Phase E — Later aerosol/PBE development (not part of first implementation)

**Status:** Not started. Intentionally deferred until the gas-phase primal and derivative path are validated.

Only begin this phase after the gas-phase primal and gas-phase derivative path is validated.

## Prompt E1 — aerosol moment transport, PBE still off

```text
Design the smallest extension to add transported aerosol moment fields while
keeping local PBE source updates disabled. Account for six physical moment
fields per mode, internal Z=Y/rho transport auxiliaries, boundary conditions,
state registration, residual fields, and transport residuals. Do not implement
microphysics, positivity clipping, vaporSink, latentSource, ADR, or objectives
in this step.
```

## Prompt E2 — PBE primal coupling only

```text
Design a primal-only extension that reproduces the standalone fractional-step
algorithm: gas PIMPLE solve, moment transport outside PIMPLE, local PBE update,
and time-lagged vaporSink/latentSource for the next time step. Identify all
non-smooth operations and explicitly state that adjoint support is not claimed.
Do not implement until standalone-versus-wrapper primal equivalence is defined.
```

## Prompt E3 — PBE derivative audit

```text
Audit the PBE/microphysics implementation for differentiability. Identify
positivity clipping, activation thresholds, freezing decisions, table lookups,
nonlinear local updates, and time-lagged coupling. Propose a staged derivative
verification sequence beginning with a reduced smooth configuration. Do not
implement ADR changes yet.
```

---

# Quick prompts for review after every Copilot patch

Use these whenever a local model proposes a change:

```text
Review this diff against the current reduced DAContrailFoam LES scope. Identify:
1. physics mismatches with CASSANDRA contrailFoam;
2. incorrect DAFoam state/residual naming;
3. any solve(), clamp(), or state modification inside calcResiduals();
4. uses of constant-property thermo that are invalid for multiComponentMixture;
5. accidental aerosol/PBE dependencies;
6. incorrect dynamicKEqn, LES, RAS, or turbulence-model assumptions;
7. generated-file modifications;
8. original/ADR/ADF build implications.
Do not edit files; report only concrete findings.
```

```text
Give the exact minimal Docker build and test commands for this patch. Build
original mode only unless ADR is explicitly requested. Include the commands to
check build status, tail the persistent log, inspect git diff --check, and avoid
simultaneous builds.
```

## Stop conditions

Stop and inspect before proceeding if any of the following occurs:

- a local model proposes `stateInfo_["volScalarStates"].append("Y")`;
- a residual calls `solve()`, clamps a state, or advances time;
- a patch uses fixed molecular weight, fixed Cp, or hand-coded Sutherland transport for the multicomponent case;
- `N2` is added as an independent state;
- code reads a `RAS` dictionary for the LES dynamicKEqn path;
- a patch substitutes kEpsilon or dummy turbulence behavior while claiming dynamicKEqn support;
- a patch claims dynamicKEqn or LES adjoint support before the dedicated turbulence residual, original-mode LES validation, and derivative checks exist;
- a patch edits `lnInclude/`, `Make/linux*`, `platforms/`, or generated case output;
- a build attempt tries ADR/ADF before explicit ADR preparation.

## Completion definition for the first milestone

The first milestone is complete only when all conditions hold:

1. `DAContrailFoam` runs the PBE-disabled, multicomponent, dynamicKEqn LES case in original mode. — **done (Stage 3, commit `6f7a7ca`; `da_long` ran to t≈0.076)**
2. The DADynamicKEqn model reproduces the standalone OpenFOAM dynamicKEqn primal behavior on a short copied case. — **done at engineering level (step gates + twin 0.06 s runs; known chaotic `nut`/`k` divergence documented; formal 0.3 s closure optional)**
3. `DAContrailFoam` reproduces standalone final and selected time-averaged fields within declared tolerances. — **done at engineering level for the t=0.06 window (`UMean` 3.5%, `TMean` 0.8%, `rhoMean` 1.0%); full-run final-field sign-off optional**
4. `U`, `p`, `T`, `rho`, `O2`, `CO2`, `H2O`, dynamicKEqn model-state fields, `nut`, and relevant thermal-turbulence fields are compared. — **done (npz cmp + patch/region statistics + `*Mean` fields; `nut`/`k` show the documented Germano amplification)
5. The state/residual system contains independent `O2`, `CO2`, and `H2O` states and derived `N2`. — **done (structure: state info, YEqn reconstruction, residual class)**
6. No aerosol/PBE, vaporSink, latentSource, or microphysics code is active in the wrapper. — **done (verified in field creation and equation includes)**
7. The implementation has a small LES regression case and a documented build/run procedure. — **pending (Phase C3); build procedure documented in `DAContrailFoam_runbook.md`**
8. ADR/ADF support is explicitly marked as future work unless separately implemented and derivative-validated. — **done (marked; AD libs stale, rebuild gated on Phase D4)**
