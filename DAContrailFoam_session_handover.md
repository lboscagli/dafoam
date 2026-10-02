# Session Handover — 2026-10-02 (Phase D4: ADR divergence root-caused and fixed; G1/G2/G3 PASS)

> Read **`DAContrailFoam_runbook.md` → "Phase D4"** for the full evidence and
> the restart commands; this file is the short-form state. The gate
> definitions are in **`DAContrailFoam_phaseD4_adr_audit.md` §9** (line ~179):
> G1 line 181, G2 line 185 (≤1e-10 residual parity), G3 line 189, G4 line 195,
> **G5 line 199 (the acceptance gate)**, G6 record. **No "ADR support" claim
> may be made until G5 passes.**

## 1. THE BLOCKER WAS FOUND AND REMOVED (root cause #6)

The ADR build linked fine but **every ADR primal diverged**, while the ADF
primal ran. The environment, the physics/assembly and the preconditioners
were all ruled out first (identity preconditioner also diverged → not DILU;
`initRes` matched at 16 digits → assembly is correct; ADR-only process with
the union env still failed → not env interposition).

Real cause: OpenFOAM's generic `gSumProd` (`Field/FieldFunctions.C:488`)
implements `sumProd` with `operator&&` — correct "dotdot" for vectors, but a
**logical-AND count** for scalars. `double`/`float` mask this because
`scalarField.{H,C}` provides `*`-based specializations for them; the codi
`scalar` (`codi::ActiveType`) had **none**, so every BiCG/BiCGStab/CG
`alpha`/`omega` became garbage → divergence of all BiCGStab solvers (the
`smoothSolver`/`GAMG` paths only use unary norms and are immune, which is
what hid it for so long).

**Fix (option B — source patch, approved):** add `sumProd` specializations
for the codi scalar in `scalarField.H` (declaration) + `scalarField.C`
(definition). The declaration is mandatory: `Field.H` includes
`FieldFunctions.C`, so consumers inline the generic and, pre-fix, there was
no symbol left for LD_PRELOAD to interpose.

**Both AD libs rebuilt and verified**, so ADR *and* ADF now carry the fix:

| lib | size | NEEDED | symbol check |
|---|---|---|---|
| `libOpenFOAMADR.so` | 24 933 632 B | `libPstreamADR.so` | 1 defined `sumProd<codi::ActiveType>`, 0 undefined |
| `libOpenFOAMADF.so` | 19 944 920 B | `libPstreamADF.so` | same |

Pre-fix backups: `/workspace/logs/backups/libOpenFOAM{ADR,ADF}.so.pre-sumProd-fix-*`.
Patched files (backups `*.orig-d4`): `scalarField.{H,C}`,
`src/OpenFOAM/Make/options` (`-lPstream` → `-lPstream$(WM_AD_MODE)` — the
image ships only the suffixed `libPstream{ADR,ADF}.so`), plus the earlier
`codedFixedValue`/`codedMixed` and functionObject `$(WM_AD_MODE)` patches.
Build logs: `logs/libOpenFOAM-ADR-relink.log`, `logs/libOpenFOAM-ADF-rebuild.log`,
`logs/libOpenFOAM-ADR-rebuild.log` (first attempt, failed at `-lPstream`).

**Build gotchas learned (both bite again):**
- `wmake` usage is `wmake [OPTION] [dir]` → **`wmake -j2 libso`**, not
  `wmake libso -j2`.
- Touching any `Make/options` regenerates the makefiles → **full 601-object
  libOpenFOAM recompile (~18 min at -j2)**, not a relink. Objects live under
  `OpenFOAM-AD/build/<WM_OPTIONS>/src/...`, **not** `platforms/`.
- `Make/files` has no AD suffix → wmake emits `libOpenFOAM.so`, which must be
  **renamed** to `libOpenFOAM{ADR,ADF}.so` (no SONAME; consumers' DT_NEEDED
  are the suffixed names).

## 2. GATES

| gate | result |
|---|---|
| **G1** ADR init + coded-BC compile + 100-step ADR primal | **PASS** — `primalFail=0`, elapsed 216.4 s, 0 singularity/NaN, fresh `dynamicCode` (9 coded BCs + `turbulenceScales_dynamicKEqn` + `fieldAverage1`), pristine DILU case. Log `logs/g1_gate_official.log` |
| **G2** residual parity original-env vs ADR-env | **PASS, bitwise identical** — 950 516 residuals, `norm2=5205369.2529319227`, `absmax=1400533.6289946034`, `max_abs_diff=0`, `elements differing: 0/950516`. Logs `logs/g2_residual_{orig,adr}.log` + `logs/g2_compare.log`, dumps `g2_residuals_{orig,adr}.npy` |
| **G3** checkpoint/replay integrity | **PASS** — 101 checkpoints (`0`…`0.00125`), `getdFScaling` nonzero exactly at indices 50..99 all `0.02=1/50` and 0 elsewhere, `primalFinalTimeIndex_` valid (Prompt-2 fatal did not fire), replayed states vary, final-state round trip **abs 5.46e-12 / rel 3.35e-16**. Log `logs/g3_replay.log`, `logs/g3_primal_final.npy` |

G2 methodology note: use `getResiduals(double*)` (`DASolver.C:1238`) — it
returns full-precision doubles. `calcPrimalResidualStatistics("print")` prints
only 6 significant digits and **cannot** be used for a ≤1e-10 gate.

G3 case changes (both required): `writeControl timeStep; writeInterval 1`
(the write sits *inside* the time loop and `writeAdjStates` is gated by
`runTime.writeTime()`, so the inherited `adjustableRunTime/0.00125` produced
exactly one time dir — no checkpoint set existed) and `writePrecision 16`
(the case ships 10; the ≤1e-12 round-trip bar is unreachable at 10 digits).

## 3. FINDING THAT BLOCKS G4/G5 — `0/` is missing `phi`

`readStateVars` `MUST_READ`s every state, and the first three steps need
time levels at/before t=0, which only exist in the hand-provided `0/` folder
(the solver never writes time 0). Observed fatal:
`cannot find file "/tmp/stage3/d4_g3/0/phi"`. DAFoam's own `mphys`
`readZeroFields` path calls `readStateVars(0.0, deltaT)` before every primal
(`mphys_dafoam.py:1487`), so **the unsteady adjoint cannot run on these
cases until `0/` holds every state field (at least `phi`)**. There is no
Python write API to fix it from the outside. G3 therefore replays the 98
checkpoints whose `t-2·dt` is also a checkpoint (t ≥ 3·dt).

Fix: `cp /tmp/stage3/d4_g3/0.00125/phi <case>/0/`, re-run G3, and require
`0/ missing state fields: []`.

## 4. WHERE EVERYTHING LIVES

- Repo: `~/dafoam-source-docker/work/src/dafoam` (git root; branch
  `feature/contrailfoam-primal-wrapper`; `origin = lboscagli/dafoam`,
  `upstream = mdolab/dafoam`).
- **Gate scripts are now committed under `d4_gate_scripts/`** (the originals
  in the host `work/scripts/` are what the container sees as
  `/workspace/scripts` — edit the host copies, then re-copy).
- Runbook: `DAContrailFoam_runbook.md` (Phase D4 section has the full
  evidence tables, root-cause #6 write-up, restart commands, gotchas).
- Audit / gate definitions: `DAContrailFoam_phaseD4_adr_audit.md` (§6 line 93
  = checkpoint mechanism, §7 = non-smooth inventory, §9 line 179 = gates).
- Cases are **container-local** (`/tmp/stage3/{d4_g1,d4_g3,d3_base}`) and are
  lost if the container is recreated; the runbook documents how to rebuild
  `d4_g1`.
- Logs: host `work/logs/` (= `/workspace/logs` in the container).

## 5. NEXT STEPS (in order)

1. **Fix `0/`** on the G4/G5 case (see §3) — everything downstream depends on it.
2. **G4** — short-window `meanTOutlet` derivative for an existing
   differentiable input, plus instrument audit §7's non-smooth inventory for
   activity in that window (known candidates: `Y.clamp_min` at
   `YEqnContrail.H:39,46`, `max(KK,0)` at `DADynamicKEqn.C:32`).
3. **G5 (acceptance)** — fixed-dt window `endTime 0.012`, objective
   `meanTOutlet`, control = jet velocity scale α, central difference with
   h ∈ {0.005, 0.01, 0.02}, `|dJ/dα|_AD − |dJ/dα|_FD| / |dJ/dα|_FD ≤ 1e-4`
   under D3's R1–R4. **Only after this may "ADR support" be claimed.**
4. **G6** — record G1–G5 evidence, flip the prompt-pack ADR/ADF bullets to
   verified.

Known-but-unreachable gaps (none affect the current cases): AD-tree
`codeStream.C:174`, `CodedFunction1.C:133`,
`codedFixedValuePointPatchField.C:119` (→libOpenFOAM) and `CodedField.C:135`
(→libmeshTools) still have unsuffixed lib references.

# Session Handover — 2026-10-01 (Phase C closed, D1–D4 done, D3 Tier-1 executed and recorded)

## 1. ARCHITECTURAL STATE

DAFoam v5.1.1 fork (branch `feature/contrailfoam-primal-wrapper`; code
commits `f33eab6` (D2 guard + record) and `fef5cc5` (Prompt 2 fail-fast
errors); D3/D4 documents committed in this changeset; git root = this
directory) wraps the CASSANDRA gas-phase LES `contrailFoam` (dynamicKEqn,
species O2/N2/CO2/H2O, 280k cells) as `DAContrailFoam` in original mode only
(`COMPILE_DAFOAM_NOAD=1`; container `dafoam-dev`; env `source
/workspace/src/cassandra/dafoam_env_doc/load_cassandra_dafoam_env.sh`).
Phases A–C closed (primal validated at engineering level); Phase D: D1 design,
D2 implementation + validation, Prompt 2 error hardening (verified by five
failure-case runs), D3 sensitivity plan, D4 ADR audit, and **D3 Tier-1
execution** all delivered (results + gates in runbook "Phase D3 Tier-1
execution"); remaining: Tier-2 window [0.1, 0.3] (deferred by user decision
2026-10-01) and ADR/ADF rebuild + FD validation on explicit request.

## 2. CURRENT CODE (ESSENTIAL SNAPSHOT ONLY)

Committed in this session (working tree clean afterwards; further commits only
on explicit request):
- `f33eab6` "Add mean objective bounds guard and phase D2 record" —
  `src/adjoint/DASolver/DASolver.C` bounds guard in `calcAllFunctions` (grow
  `functionTimeSteps_[idxI]` to `listIndex+1` zero-filled when
  `listIndex >= size`; skip store when `listIndex < 0`), plus
  `DAContrailFoam_runbook.md` (Phase D2 section: window/reset/finalization/
  restart behavior, validation numbers), `dafoam_contrailfoam_copilot_prompt_pack.md`,
  `DAContrailFoam_phaseD1_mean_objective_design.md`,
  `DAContrailFoam_session_handover.md`.
- `fef5cc5` "Fail fast on mean objective misuse" (Prompt 2) —
  `DASolver.C` (`getTimeOpFuncVal`/`getdFScaling` fatals for
  `primalFinalTimeIndex_ == 0` and unknown function name),
  `DATimeOpAverage.C`/`DATimeOpFinal.C`/`DATimeOpMax.C` empty-window guards,
  `DAFunctionPatchMean.C` `areaSum_ <= 0` fatal; runbook "Phase D2 hardening"
  section with the five verification runs.
- this changeset: `DAContrailFoam_phaseD3_sensitivity_plan.md`,
  `DAContrailFoam_phaseD4_adr_audit.md`, runbook/prompt-pack/handover status.
- Not in git: `work/scripts/stage3_primal_compare.py` (host, not a git repo) gained `--func` → injects `meanTOutlet` dict and prints `STAGE3: func meanTOutlet = ...`.

D2 validation result (runbook "Phase D2" has full record): 0.012 s copied case, 266 steps, `evalFunctions` 208.0155615012 K vs fieldAverage `TMean` outlet areaAverage 208.0159243 K → difference −3.628e-4 K = −1.744e-6 rel, **entirely dt-weighting** (DATimeOpAverage unweighted, fieldAverage dt-weighted, same sample set t1..tN); guard-trigger run (60-slot list, 120 steps) survived with exact readback. Build log: `logs/dafoam-Allmake-20261001-062527.log`.

D2 recon findings (existing machinery — no new classes needed):
- `src/adjoint/DAFunction/DAFunctionPatchMean.C`: area-weighted patch mean; dict: `type:"patchMean"`, `source:"patchToFace"`, `patches:[...]`, `scale`, `varName`/`varType`/`index` (index readEntry required even for scalar).
- `src/adjoint/DATimeOp/DATimeOpAverage.C`: unweighted step-mean, `dFScaling=1/N`.
- `DASolver::getTimeOpRange` (DASolver.C:424): window = trailing `nStepsFrac` (default 0.2) of steps → emulate [t0,t1] by `endTime=t1`, `nStepsFrac=(t1-t0)/endTime`.
- **Bug to fix first:** `functionTimeSteps_[idxI][listIndex]=functionVal` (DASolver.C:369) unguarded; list sized `round(endTime/deltaT)` (DASolver.C:584) → OOB if dt shrinks or case restarts at nonzero timeIndex. Guard: resize-on-demand to `listIndex+1` (zero-fill) before the store.
- Wiring already live: `setDAFunctionList()` via `createAdjoint.H:31`; per-step `calcAllFunctions(printToScreen_)` at `DAContrailFoam.C:178`; `primalFinalTimeIndex_` set at `DAContrailFoam.C:208`; empty `function` dict → early return (DASolver.C:342).
- Python access: `PYDAFOAM.evalFunctions` → `solver.getTimeOpFuncVal(name)` (pyDAFoam.py:926–933); `function` is a standard pyDAFoam option (line 882).
- Reference ground truth: case `system/controlDict` has `functions/fieldAverage1` (line 377: `timeStart 0.0`, `executeInterval 1`, `writeControl writeTime`, fields U/T/rho/Sw_local mean+prime2Mean, `base time`); `0.06/TMean` exists in `da_long`. OF fieldAverage IS dt-weighted (`fieldAverageItem.C:137` `totalTime_ += deltaTValue`); DAFoam `DATimeOpAverage` is NOT → quantify difference in validation (or add dt-weighting if user wants exactness). Unverified: whether fieldAverage averages boundary fields (no `boundaryField` hits in fieldAverageTemplates.C) — check `0.06/TMean` header for boundary values before relying on patch comparison.
- Driver: `work/scripts/stage3_primal_compare.py` (host `/Users/lbosc/dafoam-source-docker/work/scripts/`, not a git repo) `mode da` builds `daOptions={"solverName":"DAContrailFoam"}` → inject `"function": {...}` dict there. Container cases at `/tmp/stage3/{da,of}_long` (+ `da_s2`, `da_t*`, `of_t*` step-gate dirs). `da_long.log`: startup ~6 s, run stopped at t≈0.076; validation run must copy case and shrink `endTime` (e.g. 0.006–0.012, `writeInterval` matching).
- Build after any C++ edit: `docker exec -e COMPILE_DAFOAM_NOAD=1 dafoam-dev bash --noprofile --norc /workspace/src/cassandra/dafoam_env_doc/build_dafoam_source.sh`; authoritative check = grep `logs/dafoam-Allmake-latest.log` for `error` (status file lies). `pkill` gotcha: use `pkill -f "contrail[F]oam"`.

## 3. NEXT STEPS (DIRECT PROMPTS)

- [Prompt 1, DONE as `f33eab6`]: D2 implementation (bounds guard, rebuild, `--func`, validation vs fieldAverage1, runbook record) — see runbook "Phase D2".
- [Prompt 2, DONE as `fef5cc5`]: error handling for the D2 objective path —
  `primalFinalTimeIndex_ == 0` fatal in `getTimeOpFuncVal`/`getdFScaling`,
  unknown-function fatal in `getTimeOpFuncVal` (existing one verified in
  `getdFScaling`), empty-window guards in `DATimeOpAverage`/`Final`/`Max`,
  `areaSum_ <= 0` fatal in `DAFunctionPatchMean`; five failure-case runs
  verified (runbook "Phase D2 hardening"). Side finding: `"patches": []`
  with `patchToFace` segfaults during `PYDAFOAM` init (unisolated, follow-up).
- [Prompt 3 / D3, DONE]: plan in `DAContrailFoam_phaseD3_sensitivity_plan.md`
  (control = jet velocity scale α in `0/U`, objective = `meanTOutlet`,
  fixed dt 1e-5, Tier-1 [0.02, 0.06] s / Tier-2 [0.1, 0.3] s windows,
  R1–R4 repeatability, G1–G7 acceptance, ≈8.6 h Tier-1 matrix).
- [Prompt 4 / D4, DONE]: audit in `DAContrailFoam_phaseD4_adr_audit.md`
  (non-smooth inventory, coded-BC runtime-compile blocker, reduceIO replay +
  `getdFScaling` time-average path, build/test gates G1–G6, FD acceptance
  1e-4 on the 0.012 s `meanTOutlet` window).
- [D3 Tier-1, DONE 2026-10-01]: run matrix executed on a user-approved
  coarse mesh (14/20/60 = 78,960 cells, dt 1.25e-5, 4800 steps, window
  [0.02, 0.06] = last 3200 samples, `nStepsFrac` 0.6667; 6 matrix runs +
  series cross-check run + 0.15 s diagnostic probe). Verdict: G1/G2/G3/G7
  PASS (bit-identical determinism, plumbing 5.7e-13, Co 0.3146, all
  `primalFail=0`); G4/G5/G6 FAIL — window is non-stationary (plume arrival
  ramp + ~0.07 s outlet oscillation), so per plan the sensitivity in this
  window is reported **not identifiable at this window length**; no gradient
  claim. Secondary plateau FD [0.01, 0.03]: slopes 0.028061/0.028074 agree
  to 0.048 %, dJ/dα ≈ +2.91e-04 K per m/s. Full record: runbook "Phase D3
  Tier-1 execution"; analysis: `work/scripts/d3_sensitivity_analysis.py`.
- Next: (a) optional Tier-2 [0.1, 0.3] on explicit request (endTime 0.3,
  24000 steps, same `nStepsFrac`; expect ~4 h contended; σ_J projection
  0.028 K → G5 may still fail — re-measure, don't assume); (b) on explicit
  request, ADR/ADF rebuild starting at D4 gate G1 through G5 (note:
  `libDASolverADR.so` has 28 `DAContrailFoam` symbols but predates the
  D2/Prompt-2 edits; `libDASolverADF.so` has none — joint rebuild required).
  Optional D2 follow-ups: dt-weighted `timeOp`, recompute-from-t0 on restart.
  Do NOT commit unless the user explicitly asks.
