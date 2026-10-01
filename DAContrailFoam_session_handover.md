# Session Handover — 2026-10-01 (Phase C closed, D1+D2 done, D2 committed as f33eab6)

## 1. ARCHITECTURAL STATE

DAFoam v5.1.1 fork (branch `feature/contrailfoam-primal-wrapper`, head `f33eab6` — local, ahead of origin by 2, not pushed; git root = this directory) wraps the CASSANDRA gas-phase LES `contrailFoam` (dynamicKEqn, species O2/N2/CO2/H2O, 280k cells) as `DAContrailFoam` in original mode only (`COMPILE_DAFOAM_NOAD=1`; container `dafoam-dev`; env `source /workspace/src/cassandra/dafoam_env_doc/load_cassandra_dafoam_env.sh`). Phases A–C are closed (primal validated at engineering level: bit-exact inlet, means ≤3.5%, outlet-BC suspicion cleared, VTK suffix explained as write-metadata); Phase D: D1 design doc written, **D2 implemented, validated and committed (`f33eab6`)** — working tree clean, next is D3 (or Prompt 2 error hardening).

## 2. CURRENT CODE (ESSENTIAL SNAPSHOT ONLY)

Committed as `f33eab6` "Add mean objective bounds guard and phase D2 record" (working tree clean; further commits only on explicit request):
- `src/adjoint/DASolver/DASolver.C` — bounds guard in `calcAllFunctions` (grow `functionTimeSteps_[idxI]` to `listIndex+1` zero-filled when `listIndex >= size`; skip store when `listIndex < 0`).
- `DAContrailFoam_runbook.md` (Phase D2 section: window/reset/finalization/restart behavior, validation numbers), `dafoam_contrailfoam_copilot_prompt_pack.md`, `DAContrailFoam_phaseD1_mean_objective_design.md`, `DAContrailFoam_session_handover.md`.
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
- [Prompt 2]: "Now implement error handling for the D2 objective path: fatal-error when `getTimeOpFuncVal`/`getdFScaling` is called with `primalFinalTimeIndex_==0` (primal not run), when `functionName` is not found (already fatal — verify), and when the `function` window is empty (`iEnd < iStart` → clear error instead of div-by-zero in DATimeOpAverage); guard `DAFunctionPatchMean` against `areaSum_==0`. Verify with one short run each for the failure cases."
- [Prompt 3, later]: D3 sensitivity plan (design only, no ADR) → Prompt 4: D4 ADR audit → only on explicit user request rebuild ADR/ADF + FD-validate one short-window mean objective. Optional D2 follow-ups: dt-weighted `timeOp` for exact fieldAverage agreement; recompute-from-t0 on restart before any adjoint use. Do NOT commit unless the user explicitly asks.
