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
- **Primal equation set (structurally present, not yet validated):**
  `solvePrimal()` assembles the momentum equation
  (ddt + multivariate div + `muEff` laplacian + fvOptions) and includes the
  enthalpy, species, and pressure-correction equations ported from CASSANDRA
  `hEqn_aero.H`/`YEqn_aero.H`/`pEqn`, with multivariate `mvConvection`,
  `muEff/Sc` species diffusion, and inert-species reconstruction
  `Y_N2 = 1 - ΣY_k`. The PIMPLE loop structure runs, but it is still the
  single-shot skeleton: no outer time loop, equation order and auxiliary-field
  updates differ from CASSANDRA, and residual evaluation currently aborts on a
  dimensions mismatch — i.e. it compiles, wires up, and marches code, but is
  **not yet a physically validated transient primal** (Stage 3).

In short: discovery, construction, field creation, turbulence selection, and the
entire adjoint plumbing are complete and verified in original mode; the physics
time loop and residual/adjoint numerics are the remaining work.

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

The codebase is now in a small but meaningful incremental state:

- DAFoam can recognize the new solver family
- state metadata exists for the no-PBE gas-phase path
- the residual class contains the reduced compressible species transport structure
- aerosol/PBE support is still intentionally deferred

## Recommended next step

The next reviewable patch should continue the reduced multicomponent gas-phase work without broadening scope:

1. tighten the exact active species list handling in the state metadata
2. verify the residual equation matches the exact standalone contrailFoam `YEqn_aero.H`
3. validate a small no-PBE case with `aerosolActive false`
4. only then consider turbulence-model support or subsequent aerosol extensions

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

## Known remaining items (not part of this patch)

- Stage 3: `solvePrimal` still runs a single skeleton PIMPLE loop; the outer
  `runTime` loop, equation ordering vs CASSANDRA (`UEqn`→`YEqn`→`hEqn`),
  `rhoEqn`, `storePrevIter`, `pressureControl.limit`, CourantNo/`setDeltaT`,
  turbulence `correct()` and field writing are not yet parity-checked.
- Residual evaluation probe (`calcPrimalResidualStatistics`) aborts with a
  dimensions mismatch (`[1 2 -3]` vs `[1 -1 -3]`) inside residual math — Stage 3 scope.
- `DAField.C:1137` hard-codes `turbDict.subDict("RAS")`, but only fires when
  `primalBC` contains `useWallFunction`; default `primalBC = {}` keeps it dormant.
- ADR/ADF builds must be rebuilt before reverse/forward AD runs can use
  `DAContrailFoam` (original-mode `COMPILE_DAFOAM_NOAD=1` skips them).

## Notes

The implementation should stay intentionally narrow and reviewable. The first milestone is not a full DAFoam contrail solver; it is the reduced gas-phase path that matches the validated CASSANDRA no-PBE physics.
