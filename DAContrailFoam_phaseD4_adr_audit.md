# DAContrailFoam — Phase D4: ADR enablement audit (checklist + build/test plan)

Status: audit complete (documentation only — **no code edits, no AD rebuild**
in this phase, per the prompt-pack gate: the ADR/ADF rebuild happens only on
explicit request, after this checklist is agreed).

**Non-claim:** no ADR/ADF support is claimed for `DAContrailFoam` anywhere in
this repository until gate **G5** (finite-difference validation of one
short-window mean objective) passes. The primal LES evidence from Phase C is
explicitly *not* derivative evidence.

## 1. Files audited

- `src/adjoint/DASolver/DAContrailFoam/*` (`DAContrailFoam.C`, `UEqn/H/Y/p/
  rhoEqnContrail.H`, `createFieldsContrail.H`, `createRefsContrail.H`)
- `src/adjoint/DAResidual/DAResidualContrailFoam.C`
- `src/adjoint/DAStateInfo/DAStateInfoContrailFoam.H`
- `src/adjoint/DAModel/DATurbulenceModel/DADynamicKEqn.{H,C}`
- `src/adjoint/DAFunction/*`, `src/adjoint/DATimeOp/*` (the objective path,
  incl. the Phase D2/Prompt-2 guards)
- multicomponent thermo path: `createFieldsContrail.H` →
  `psiReactionThermo::New`, `multiComponentMixture` species
  (`O2 N2 CO2 H2O`), `constant/thermophysicalProperties`
- external libs: OpenFOAM-AD (`dynamicKEqn`, `fvOptions`, `fvPatchField`
  coded-BC machinery, `topoSetSource`), PETSc, the runtime `dynamicCode`
  compiler, `libDASolverADR.so` / `libDASolverADF.so`
- driver/access path: `dafoam/pyDAFoam.py`, `dafoam/mphys/mphys_dafoam.py`,
  `src/adjoint/DAInput/*`

## 2. Build mechanics (verified from `Allmake`)

- `Allmake:27-41`: `COMPILE_DAFOAM_NOAD` **unset** → ADR mode
  (`WM_AD_MODE=ADR` in `OpenFOAM-AD/etc/bashrc`), then `COMPILE_DAFOAM_ADF`
  set additionally → ADF mode. Our usual rebuilds set
  `COMPILE_DAFOAM_NOAD=1`, so only the original lib is refreshed.
- Current libs (container): `libDASolver.so` **Oct 1 07:38** (current head),
  `libDASolverADR.so` **Sep 21**, `libDASolverADF.so` **Sep 15** → both AD
  libs predate Stage 1+2, Stage 3, Phase A–D2 and contain **no
  `DAContrailFoam`** in the ADF case (prompt-pack finding) → any ADR/ADF run
  today is invalid for this solver.
- The AD build recompiles the *same* DAFoam sources against OpenFOAM-AD
  (codi-transcribed `scalar`), so every source-level incompatibility below is
  a build-time or runtime error in that pass.

## 3. Component-by-component audit

| component | AD-ready? | findings / required original→ADR change |
|---|---|---|
| `DAContrailFoam::solvePrimal` time loop | mostly | `setDeltaT.H` (:128) is dt control — must be identical between primal and replay (fixed dt, D3 §4); `runTime.write/writeNow` are I/O only; `validateStates` early `return 1` (:166-172) **aborts the horizon** — must be guaranteed not to fire in ADR/FD runs (it silently changes dJ's window); `primalFinalTimeIndex_` is only set on a completed loop (:208) → replay must set/inherit it (see §6). |
| `UEqnContrail.H` | yes | no clips/branches found beyond `pimple.momentumPredictor()` (config-constant branch). |
| `YEqnContrail.H` | **non-smooth** | `Y[speciesI].clamp_min(0)` (:39) and `Y[inertIndex].clamp_min(0)` (:46): kink at Y=0. Active every step if any species approaches 0 (O2/CO2/H2O are near-zero in places) → subgradient in AD, possible FD mismatch. Verify activity (G4b) before the FD gate. |
| `pEqnContrail.H` | **non-smooth** | `pressureControlPtr_->limit(p)` (:71): clamp + `if` branch. `pressureControl` in this case must be checked for `pRefNeeded`; if inactive the branch is constant-off (harmless), if active the kink is on the differentiation path. |
| `hEqnContrail.H` / `rhoEqnContrail.H` | yes | energy/continuity forms, no clipping found. |
| `DAResidualContrailFoam.C` | **non-smooth** | `Y_.clamp_min(0)` (:371, :404) in the residual/BC path (same kink as above); `correctBoundaryConditions()` (:385-405) calls `U_,p_,T_,O2_,CO2_,H2O_,Y_` → **this is the path that evaluates the inlet `codedFixedValue` BCs inside residual assembly** (see §4). |
| `DAStateInfoContrailFoam` | n/a | state bookkeeping only (`9*nCells + nFaces` states, `N2` reconstructed). No arithmetic. |
| `DADynamicKEqn.C` | **non-smooth, continuous** | inventory: `max(KK, 0)` in `MM` (:32); soft-clamps `0.5*(mag(Ck)+Ck)` (:47) and `0.5*(mag(Ce)+Ce)` (:64); `KK.clamp_min(SMALL)` (:76, :258); `sqrt(k_)` in `correctNut` (:86) and `kEqn` (:270) (safe only while k>0 — the clamp above enforces it, but that clamp is itself the kink); division regularized by `VSMALL` (:41-42) → smooth. **Good news:** the dynamic procedure is pointwise Germano with filtered ratios — no `LUsolve`/global least-squares/rank-dependent solve anywhere in the file, so the classic dynamic-model AD blocker (rank-deficient least squares) is absent. `filter_/simpleFilter_` are linear. |
| multicomponent thermo | mostly yes | `thermoType`: `perfectGas` + `specie` + `sutherland` transport (no polynomial/spline tables, no `chemistry` directory, no reactions) → Sutherland law and perfect-gas EOS are smooth for T>0. Watch for `mag(T)`-style guards upstream in OpenFOAM thermo (verify at build). |
| objective path (`DAFunctionPatchMean`, `DATimeOp*`, `DASolver` time-op code) | yes | pure arithmetic + label indexing; D2/Prompt-2 guards branch on label comparisons and on `areaSum_ <= 0.0` (a value test on geometry-only data, executed for safety, not differentiated). `getTimeOpRange` uses `round(nStepsFrac * scalar(...))` on a value cast from the function dict — under ADR `scalar` may be a codi type here; if so `round(ad)` will not compile → **anticipated build fix: cast the window arithmetic to `label`/`double` explicitly** (window indices are never differentiated). |
| coded BCs (`0/U`, `0/T`, `0/H2O`, `0/*` outlet `inletOutlet`) | **top risk** | `codedFixedValue` is compiled at **runtime** by `CodeStream` (case has `dynamicCode/`) using `::sqrt`, `::atan2`, `::tanh`, string-literal `if` branches. Under ADR these bodies are compiled during case setup — they must (a) find the AD headers so `scalar` resolves to the codi type, (b) compile against AD-safe math overloads, and (c) link. **Unverified — no evidence either way yet.** Because `DAResidualContrailFoam::correctBoundaryConditions()` runs them inside residual assembly, a failure here blocks everything. Mitigation order: 1) try the ADR runtime compile as-is; 2) replace the profiles with statically compiled AD-capable BCs (or pre-built `fixedValue` profile fields); 3) expose the control through `DAInputPatchVar`/`patchVariable` input instead of editing BC source (§5). |
| `pyDAFoam` / mphys plumbing | yes, conditional | `evalFunctions` returns `getGradient()` under `CODI_ADF` and `getValue()` under `CODI_ADR` (`DASolver.C:482-488` already handles both); `useAD mode ∈ {forward, reverse}`; inputs available for a BC control: `DAInputPatchVar`, `DAInputPatchField`, `DAInputPatchVelocity` (registered input types) — the FD file-edit control of D3 would be replaced by one of these for the AD run. |
| `validateStates` / field validators | primal-only | comparison-based abort checks; never differentiate through them, but they can truncate a perturbed run differently (FD hazard). Add to the run gate: `primalFail=0` for *every* ± run (already in D3 G7). |

## 4. The coded-BC problem (most likely blocker)

Chain of evidence:

1. Inlet/outlet conditions are `codedFixedValue` / `inletOutlet` with
   runtime-compiled `code` blocks (`0/U`, `0/T`, `0/H2O`).
2. `DAResidualContrailFoam::correctBoundaryConditions()`
   (`:385-405`) calls `correctBoundaryConditions()` on every state during
   residual evaluation — i.e. the BC bodies execute on the differentiation
   path, not just in the primal.
3. `dynamicCode/` in the case directory confirms the CodeStream
   compilation happens per case (and on every BC edit, as D3 exploits).

Required verification (first item of build gate G1): does an ADR-mode
`initSolver()` compile and load the dynamic BC code? Outcomes:
- compile error → mitigation 2/3 above;
- loads but compiled against non-AD `scalar` → **silent derivative loss**
  (worst case; must be excluded by a residual-parity test G2, not assumed);
- works → record the flags that made it work.

## 5. Control exposure for the AD run (decide before G5)

D3's control is a file edit (`0/U` constant) — fine for FD, invisible to the
adjoint. For the AD side, the same α must become a registered input:
`DAInputPatchVar`/`patchVariable`-style input mapped to the inlet `U`
profile amplitude, or (simpler, still valid for validation) keep the FD
control as-is and differentiate a *different* available parameter. The
decision affects only the FD gate, not the build: pick **patch-variable
input for α** so AD and FD perturb the same physical quantity by the same h.

## 6. Unsteady checkpoint/replay and the time-averaged objective

Mechanism actually present in this codebase (verified):

- **Checkpoints:** `reduceIO` (default `True`, `unsteadyAdjoint.reduceIO`)
  makes `solvePrimal` call `writeAdjStates(reduceIOWriteMesh_, …)` every step
  (`DAContrailFoam.C:189-193`), i.e. per-step state directories under the
  case — the replay source. Known omission from Phase C: reduceIO writes
  leave out `N2` (reconstructed — fine), `nut`, `alphat`, `rho`, and
  `uniform/time`'s `index`.
- **Replay:** `mphys_dafoam.py:1685-1710` loops time instances, calls
  `DASolver.readStateVars(timeVal, deltaT)` (so it replays by **time value**,
  not by the `uniform/time` index — the missing `index` entry is therefore
  expected to be harmless, but this must be *confirmed*, not assumed),
  then accumulates the objective derivative per step.
- **Time-average derivative:** for each step n it calls
  `getdFScaling(firstFunctionName, n - 1)` and adds
  `dFScaling * dFdW_n` into the accumulated `dJ` — exactly the D1 design
  (`DATimeOpAverage::dFScaling` returns 1/N; window membership from
  `getTimeOpRange`). Notes/requirements:
  - `getdFScaling` reads only the **window indices**, not the stored
    samples, so replay does not need `functionTimeSteps_` values — but it
    *does* need `primalFinalTimeIndex_ != 0`; after Prompt 2 a replay
    process without a completed primal now fails loudly (correct behavior —
    previously it returned 0 silently).
  - all functions are assumed to share one `timeOp`
    (`mphys` comment at :1700) — satisfied today (single mean objective).
  - `unsteadyAdjoint.mode` defaults to `"None"`; the replay path requires
    `"timeAccurate"` (or `"hybrid"`) to be set explicitly — currently unset
    in our driver options → **configuration gap, not a code gap**.
  - if a dt-weighted timeOp is ever adopted (D2 follow-up),
    `dFScaling` becomes `dt_n/Σdt` — a design change that must be made in
    `DATimeOpAverage` *and* kept consistent with `D3`'s fixed-dt window;
    until then fixed dt makes unweighted and dt-weighted identical.
  - checkpoint dt sequence must match the primal exactly (D3 §4 fixed dt)
    or the replayed states correspond to different times.
- **dynamicKEqn replay needs:** `simpleFilter_`/`filter_` are reconstructed
  from `U_` per evaluation (no extra stored state), `nut_` is
  `correctNut()`-regenerated; `oldTime` fields required by `ddt` schemes are
  handled by `readStateVars`/`additionalOldTime` options — verify on the
  first replay (G3).

## 7. Non-smooth operation inventory (complete list found)

| # | operation | location | active in nominal run? | AD behavior | FD risk | action |
|---|---|---|---|---|---|---|
| 1 | `Y.clamp_min(0)` | `YEqnContrail.H:39,46` | likely (species near 0) | subgradient 0/1 | moderate | instrument in G4b; if active, expect kink-level FD disagreement |
| 2 | `Y_.clamp_min(0)` | `DAResidualContrailFoam.C:371,404` | same | same | moderate | same |
| 3 | `pressureControl.limit(p)` | `pEqnContrail.H:71` | check `pRefNeeded` | branch + kink | low if inactive | verify inactive; document |
| 4 | `max(KK, 0)` | `DADynamicKEqn.C:32` | yes (near-zero TKE regions) | subgradient | low-moderate | accept; tolerance in G5 |
| 5 | `0.5*(mag(Ck)+Ck)`, `0.5*(mag(Ce)+Ce)` | `DADynamicKEqn.C:47,64` | yes | non-smooth at 0 (continuous) | low | accept |
| 6 | `KK.clamp_min(SMALL)` | `DADynamicKEqn.C:76,258` | yes | kink | low | accept |
| 7 | `sqrt(k_)` | `DADynamicKEqn.C:86,270` | yes | smooth for k>0 (guaranteed by #6) | low | guard dependency #6 |
| 8 | coded BC runtime compile + `if ("Uz_only"=="Uz_only")` string branch | `0/U` etc. | yes | branch is compile-time constant (folded) | — for branch; **high** for compile/link (§4) | G1/G2 test |
| 9 | `validateStates` abort | `DASolver.C:3868+`, `DAContrailFoam.C:166` | only on bad runs | changes horizon | high *if it fires* | require `primalFail=0` in every FD/AD run |
| 10 | `setDeltaT` / adaptive dt | `DAContrailFoam.C:128` | if `adjustTimeStep yes` | changes discretization | **high** (dt-path noise) | fixed dt (D3 §4) |
| 11 | `round()` on possibly-ad `scalar` | `DASolver.C:457` (window) | yes | **may not compile** in ADR | — | anticipated build fix (§3) |
| 12 | `if (pressureControlPtr_->limit(p))` return-bool branch | same as #3 | see #3 | branch | see #3 | see #3 |

Explicitly checked and **absent**: global least-squares/rank-dependent solves
in the dynamic model (`LUsolve`/`svd`/`inv` — none), reaction-rate tables
(no chemistry), spline/polynomial property tables (`perfectGas` +
`sutherland`), `sign()`/`round()` in the equation headers, PBE/microphysics
coupling (`vaporSink`/`latentSource` absent from the gas solver).

## 8. Original→ADR build-change checklist (expected work)

- [ ] Rebuild: unset `COMPILE_DAFOAM_NOAD`, set `COMPILE_DAFOAM_ADF=1`
      (or ADR only); full recompile against OpenFOAM-AD; archive the log as
      the evidence file (same naming convention as the original builds).
- [ ] Fix anticipated compile items: window `round()` on codi scalars
      (`DASolver.C:457`), any `double`/`word`/printf mismatches surfaced in
      `DAContrailFoam/*`, `DADynamicKEqn`, `DAFunction`, `DATimeOp`;
      codi-compatibility of `FatalErrorIn` streams in guards (labels only —
      expected fine).
- [ ] Runtime-compiled coded BCs under ADR (§4): make the dynamic code
      compile, or replace with statically compiled equivalents.
- [ ] Confirm `libDASolverADR.so` now exports `DAContrailFoam`
      registration (today's ADR lib predates it; ADF lib lacks it).
- [ ] Register the AD input for α (§5) if not already covered by an
      existing `patchVariable` input.
- [ ] Set `unsteadyAdjoint.mode = timeAccurate` (and `reduceIO` as needed)
      in the FD/AD driver options.
- [ ] Decide dt-weighted vs unweighted timeOp (D2 follow-up) *before*
      locking `getdFScaling` semantics for the adjoint.

## 9. Build/test plan (gates, in order; G5 is the acceptance gate)

- **G1 — ADR compile + init:** clean ADR build (0 errors, log archived);
  then `initSolver()` + a few primal steps on the reference case must run in
  ADR mode, including the coded-BC dynamic compile (§4). *Fail → apply §8
  mitigations and repeat.*
- **G2 — residual parity:** `calcPrimalResidualStatistics` in ADR mode vs
  the original-mode numbers on the same case/states: agreement to
  ≤1e-10 (structural equality expected; catches silent non-differentiable
  sub-paths such as a BC compiled against the wrong `scalar`).
- **G3 — replay integrity:** one short `timeAccurate` replay over the
  checkpoint set (e.g. the first 100 steps): replayed states match the
  primal states (≤1e-12), `getdFScaling` returns 1/N in-window and 0
  out-of-window, and `primalFinalTimeIndex_` is valid in the replay process
  (Prompt-2 fatals must *not* fire — that would prove the plumbing is
  wrong).
- **G4 — derivative plumbing without chaos:** (a) short-window
  `meanTOutlet` derivative for an *existing* differentiable input (e.g. a
  patch-variable input) on a fixed-dt window; (b) instrument the non-smooth
  inventory (§7) for activity in that window.
- **G5 — FD validation of one short-window mean objective (acceptance):**
  fixed-dt window `endTime 0.012`, objective `meanTOutlet` (D2-validated),
  control = jet velocity scale α (D3 definition) exposed per §5;
  central difference with h ∈ {0.005, 0.01, 0.02};
  **acceptance: |dJ/dα|_AD − |dJ/dα|_FD| / |dJ/dα|_FD ≤ 1e-4** for the h
  that satisfies D3's SNR rule, with the D3 repeatability gates R1–R4
  satisfied (fixed dt ⇒ FD and AD see the same discretization). If a §7
  clip is active in the window, either (i) re-target a window where it is
  inactive, or (ii) declare the relaxed tolerance *before* looking at the
  numbers. Only a passed G5 permits the phrase "ADR support" in this repo.
- **G6 — record:** add the G1–G5 evidence (log paths, numbers, tolerances)
  to the runbook and flip the prompt-pack ADR/ADF bullets from "stale, gated"
  to the verified state.

## 10. Open questions carried into the build

1. Do runtime-compiled `codedFixedValue` BCs link and differentiate under
   ADR? (§4 — decides whether the inlet profiles must be re-implemented.)
2. Is `pressureControl` engaged in this case (`pRefNeeded`)? (#3)
3. Do the `Y.clamp_min` clips actually activate in the validation window?
   (#1/#2 — decides the achievable FD tolerance.)
4. Does replay tolerate the reduceIO omission of `uniform/time` `index`?
   (§6 — `readStateVars(timeVal, …)` suggests yes.)
5. Which input type maps α most directly: `patchVariable` on inlet `U`, or
   `patchVelocity`? (§5)
