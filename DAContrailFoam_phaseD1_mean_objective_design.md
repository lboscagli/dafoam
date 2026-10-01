# Phase D1 — Time-averaged LES mean-observation objective (design only)

Status: design delivered per Prompt D1 (no code written). This is the first
work item of Phase D; D2 implements exactly this design.

## 1. Observable

First objective: **post-transient time-mean temperature on a fixed region**,
with velocity mean as a drop-in variant (same machinery, different field):

```
J = 1/(t1 - t0) * ∫_{t0}^{t1} [ 1/V_R * ∫_R T(x, t) dV ] dt
```

- Region R (recommended first choice): the **outlet patch** (2,800 faces,
  fixed, always written, and already the focus of the Phase C comparison).
  Area-weighted face mean. Alternative regions (cellZone, sampling plane,
  probe ball) use the same interface, only a different mask/weight builder.
- Window [t0, t1]: post-transient, e.g. t0 = 0.1 s, t1 = 0.3 s on the 0.3 s
  reference run. First validation pass uses [0, 0.06] so it can be checked
  against the already-produced `fieldAverage` `TMean` at the common write.
- Governing equations and the adjoint target remain the instantaneous
  unsteady system. Only the observation operator is averaged. Nothing in
  `UEqn`/`YEqn`/`hEqn`/`pEqn`/`rhoEqn` or the turbulence model changes.

## 2. Discrete accumulation

Per time step t_n inside the window (step already completed, after
`rho = thermo.rho()`, at the existing `calcAllFunctions` call site in
`DAContrailFoam.C::solvePrimal`):

```
dt_n   = runTime.deltaTValue()
w_n    = 1            if t0 <= t_n <= t1 else 0
S     += w_n * dt_n * meanRegion(T)      # region mean already volume/area-normalized
W     += w_n * dt_n
J      = S / W                          # W > 0 once any in-window step exists
```

- Rule: left-rectangle on the step end states (deterministic; the same states
  the adjoint sees). Midpoint variant is a one-line change if bias matters.
- Steps before t0: no contribution (transient discarded). Steps after t1: no
  contribution (window closed; J frozen).
- `region mean` for a patch: `gSum(phi-free)`: use
  `gSum(area * T_patch) / gSum(area)`; for a cellZone: `gSum(V * T) / gSum(V)`.
- Normalization is by accumulated window time W, not by step count — required
  because `adjustTimeStep` makes dt variable.

## 3. Reset, finalization, restart

- **Reset:** S = W = 0 at solver start (or when t < t0). Guarded by an
  `windowActive` flag so pre-window steps cost one comparison only.
- **Finalization:** at the first step with t_n > t1 (or at endTime if
  t1 >= endTime), J is complete and must be reported exactly once; further
  steps do not modify S/W/J.
- **Restart:** S and W are pure functions of the trajectory. Preferred design:
  *recompute, do not serialize* — on restart from a checkpoint at t_r,
  re-run the accumulation from t0 to t_r only if t_r < t1 AND the pre-t_r
  trajectory states are replayed (unsteady adjoint replay does exactly this
  walk). This keeps the checkpoint format unchanged (no new fields in
  `uniform/`). If replay-based recompute is not yet available, v1 may
  serialize `{S, W, t0, t1}` into a small scalar file written alongside
  `writeAdjStates` and read on restart — acceptable for primal-only D2
  validation, must be revisited before D3/D4 adjoint use.
- The accumulator is a **plain scalar**, not a registered field: it must never
  enter the state vector, residuals, or mesh writes.

## 4. Adjunction (why this form works for the unsteady adjoint)

Per-step contribution is linear in the observed field with a fixed
coefficient:

```
c_n = w_n * dt_n / W          (constant per step, frozen once W is frozen)
J   = Σ_n c_n * meanRegion(T_n)
```

The discrete adjoint accumulates `dJ/d state_n = c_n * d meanRegion/d state_n`
at each in-window step and then runs the usual checkpoint/reverse sweep —
i.e. the existing DAFoam unsteady adjoint loop with a per-step objective
weight. Out-of-window steps contribute zero weight. No new adjoint
equations; D4 only has to confirm that the DAFoam unsteady
checkpoint/replay machinery wraps `solvePrimal`'s loop for DAContrailFoam.

## 5. Validation gate (for D2)

- Recompute `[0, 0.06]` window means on both `da_long` and `of_long` from the
  written per-step states is not possible (only write-time states exist), so
  D2 validation is: run both solvers with `writeInterval = dt_write` small
  enough (or use the `fieldAverage` object already in the case) and compare
  `J` against the case's `TMean` at the same instant. Acceptance: relative
  difference within the sampling/statistical tolerance established by the
  twin-run spread measured in Phase C (`TMean` maxRel 0.8% over 0→0.06;
  `UMean` 3.5%), not 1e-6.
- Also assert exactness properties of the implementation itself: J after
  window equals the direct average of the collected per-step region means;
  steps outside window have zero influence; dt-weighting correct under
  variable dt (synthetic test with unequal dt).

## 6. Why mean-prediction vs mean-measurement (not averaged instantaneous error)

Averaging instantaneous squared error gives

```
E[(T - T_meas)^2] = (mean bias)^2 + var_T + var_meas - 2 cov
```

- The variance terms do **not** shrink as the window grows; J would be
  dominated by irreducible turbulent fluctuation that no deterministic mean
  flow can match, so J is large even when the mean prediction is excellent.
- Its gradient is dominated by variance-reduction directions (eddy
  amplitude/phase), which are chaotic-noise-dominated: finite-difference and
  adjoint estimates would have disastrous signal-to-noise.
- Comparing mean-to-mean isolates the systematic observable that the mean
  flow actually controls, converges as O(sigma/sqrt(T_window/tau_int)), and
  gives a differentiable, FD-verifiable target for D3/D4.
- The measurement side must be averaged the same way (mean of the measured
  series over an equivalent window), so both sides estimate the same
  expectation.

## 7. Planned interface (D2 scope, listed for review before coding)

- Objective name: `meanTOutlet` (registry key in `daOption` `objective`),
  parameters: `field (T|U_mag)`, `regionType (patch|cellZone|box)`,
  `regionName`, `t0`, `t1`, `weighting (area|volume|none)`.
- Touch points: DAFoam objective-function list construction (already created
  by `createAdjoint.H`), `DAContrailFoam::solvePrimal` `calcAllFunctions`
  call (already present), a small accumulator struct local to the objective,
  reporting through the existing function print/write path.
- No changes to `DAResidual*`, `DAStateInfo*`, equations, mesh, or case
  setup beyond optionally enabling `fieldAverage` for cross-check.

## 8. D2 recon addendum (2026-09-30, mapping design → existing machinery)

The machinery above already exists in-tree; D2 is configuration + one guard,
not a new class:

- Region mean: `DAFunctionPatchMean` (`type patchMean`, `source patchToFace`,
  `patches`, `scale`, `varName/varType/index`) — area-weighted, per-step.
- Accumulation/window: `DATimeOpAverage` + `DASolver::getTimeOpRange`
  (`DASolver.C:424`) — **unweighted** step-mean over a **trailing**
  `nStepsFrac` window (default 0.2). Absolute [t0,t1] emulated by running
  `endTime = t1` with `nStepsFrac = (t1-t0)/endTime` (nStepsFrac=1.0 when
  t0=0). Reset/finalization come for free: `functionTimeSteps_` is sized once
  per run (`round(endTime/deltaT)`, DASolver.C:584) and `getTimeOpFuncVal`
  reads only `[startIdx, endIdx]` with `endIdx = primalFinalTimeIndex_-1`.
- Gaps vs this design: (a) dt-weighting absent — OF `fieldAverage` IS
  dt-weighted (`fieldAverageItem.C:137`), so validation must quantify the
  unweighted-vs-weighted difference (or a weighted timeOp is added if the
  user wants exactness); (b) no bounds guard at `DASolver.C:369`
  (resize-on-demand needed for adaptive dt / restarted timeIndex);
  (c) restart recompute-from-t0 is not implemented (per-step values are
  in-memory only) — revisit before D3/D4 adjoint use.
- Wiring already live: `setDAFunctionList()` (`createAdjoint.H:31`), per-step
  `calcAllFunctions` (`DAContrailFoam.C:178`), `primalFinalTimeIndex_` saved
  (`DAContrailFoam.C:208`), python `evalFunctions` → `getTimeOpFuncVal`.
- Ground truth in-case: `system/controlDict` `functions/fieldAverage1`
  (`timeStart 0`, `executeInterval 1`, mean+prime2Mean, `base time`).
  Open question for validation: whether fieldAverage boundary values are
  averaged (needed for patch comparison) — verify `TMean` boundary data
  first; else compare via `postProcess -func surfaceFieldValue` (areaAverage)
  on `T` per step or fall back to an internal-field region comparison.
