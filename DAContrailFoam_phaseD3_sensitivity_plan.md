# DAContrailFoam — Phase D3: LES sensitivity & identifiability plan
# (one control, one objective; no ADR; plan document — no solver code)

Status: plan complete (user waived the approval gate). Execution of the run
matrix below is a separate, explicitly costed step.

## 1. Scope and non-goals

- **In scope:** a single scalar control, a single time-averaged objective,
  finite-difference (FD) sensitivity protocol, repeatability checks,
  turbulent sampling-uncertainty budget, acceptance criteria, run matrix.
- **Out of scope:** ADR/ADF (Phase D4 gate), multi-control identifiability,
  aerosol/PBE controls (Phase E), any solver/driver code changes for this
  phase (the FD control is a one-number edit in `0/U`, run with the existing
  `stage3_primal_compare.py --func`).

## 2. Control: inlet jet velocity scale α (chosen)

**Definition.** α multiplies the jet core velocity constant in the
`0/U` inlet `codedFixedValue` BC (`UProfile`, code block line
`const scalar jetVal = 96.4;`): perturbed value `96.4*(1 ± h)`.
Ambient velocity (`ambVal = 1.0`) is unchanged. Physical meaning:
α = 1 is the reference jet; α = 1+h is a +h% jet velocity.

- Runs are restarted from the same `0/` directory each time; the coded BC is
  recompiled at runtime (`dynamicCode/`), so **no solver rebuild is needed**
  per perturbation — only a `sed` edit of the constant in `0/U`.
- The BC also contains a deterministic angular perturbation
  (`sin(8θ + 4η)`, amplitude 1%) that is scaled by `UzBase`, so it scales
  consistently with α.

**Why not ambient H2O mass fraction (the alternative in the D3 prompt).**
Verified facts: `0/H2O` inlet has `jetVal = 0.0; ambVal = 0.0`,
`internalField uniform 0.0`; `constant/` has **no chemistry directory** and
the gas solver has **no microphysics coupling** (`vaporSink`/`latentSource`
are absent from `DAContrailFoam/*.H`). H2O is therefore a passive scalar with
zero background concentration in the current gas-phase-only configuration:
the expected response of T to a ppm-level humidity change is below the
sampling noise (§7), i.e. the control would be *practically unidentifiable*
here for structural reasons, not statistical ones. The H2O control becomes
viable after Phase E couples the vapor sink; keep it on the deferred list.
A dry-run confirmation (one ±1% H2O pair) can be added cheaply after the
velocity study if a null-result baseline is wanted.

## 3. Objective: `meanTOutlet` (validated in D2)

J(α) = area-weighted mean of `T` over the `outlet` patch (2800 faces,
0.2824527133 m²), averaged in time over a post-transient window [t0, t1]:

- per-step sample: `DAFunctionPatchMean` (`type patchMean`,
  `source patchToFace`, `patches ["outlet"]`, `varName T`) — path validated
  against `fieldAverage1` in Phase D2 (agreement to dt-weighting, ≤1.74e-6
  relative);
- time average: `timeOp average` with `nStepsFrac = (t1-t0)/t1` and
  `endTime = t1`, which reproduces exactly the trailing window [t0, t1] when
  dt is fixed (§4) — the step-fraction vs time-fraction caveat in D1/D2
  disappears because dt is constant;
- extraction: `stage3_primal_compare.py --mode da --case … --func` prints
  `STAGE3: func meanTOutlet = …`; per-step values and `Time =` lines in the
  log give an independent recomputation of both the unweighted and (with
  fixed dt identical) dt-weighted mean.

## 4. Run policy (identical for every run in the matrix)

| item | setting | reason |
|---|---|---|
| time stepping | `adjustTimeStep no; deltaT 1e-5` (pilot-gated, see below) | with `adjustTimeStep yes` (`maxCo 0.4`, observed dt 1.2e-5…5.15e-5) the dt sequence becomes state-dependent: ±α runs take different dt paths and the FD difference mixes sensitivity with time-discretization noise. Fixed dt makes the sample grid identical by construction (R2). 1e-5 divides every window boundary in §5 into whole steps, and is ≤ the smallest adaptive dt observed (1.2e-5; the nominal first step already ran at 1e-5), so it is Co-safe by construction. |
| pilot gate | one run at the chosen fixed dt must reach t1 with `primalFail=0` and Co ≤ 0.4 (parse CourantNo from the log) | 1.2e-5 is the minimum adaptive dt observed in D2, so it is Co-safe; a larger fixed dt is allowed only if the pilot shows margin — it cuts wall time proportionally. |
| parallelism / binary | serial, single binary (commit `fef5cc5` or later) | reduction order must be identical across runs (determinism check R1). |
| I/O | `writeInterval = t1` (single write), `reduceIO` default | one checkpoint set per run; `fieldAverage1` left enabled for diagnostics (cheap, does not affect the primal). |
| output | run log + `--out` npz + `0/U` copy of the perturbed constant | provenance for every run. |

## 5. Timescales and windows

Measured (case geometry/BCs from `blockMeshDict`, `0/*`; D2 run logs):

- domain: 1.4 m long, inner square 0.054 m, O-grid to r = 0.3 m; mesh 280k
  cells; jet radius 9 mm, top-hat core 96.4 m/s, ambient 1.0 m/s,
  T_jet 322 K, T_amb = internalField 208 K;
- **jet transit time** L/U = 1.4/96.4 ≈ **14.5 ms** (ambient would need
  1.4 s; internal field already equals the ambient state, so only the jet
  plume is a transient);
- D2 window (0–0.012 s, 266 samples): outlet T departs 208 K at
  t ≈ 3.9 ms (startup adjustment), mean 208.015561 K, σ_inst = 1.28e-2 K
  (dominated by the arrival ramp), tail-half σ = 3.1e-5 K, lag-1 acf 0.999,
  **τ_int ≈ 0.87 ms**, σ_J(0–0.012 s) = 7.9e-4 K;
- Phase C: FP-level seeds grow to k-rel-diff 1.35e-2 within 3 steps; twin
  runs diverge monotonically to 0.076 s → the flow is strongly chaotic and
  sampling error, not roundoff, dominates long-window means.

**Windows (two tiers):**

- **Tier 1 (screening): t0 = 0.020 s, t1 = 0.060 s** — starts after
  >1 transit and after the startup adjustment, ends at the Phase C common
  write time (comparable flow state already analysed). Window T_win = 40 ms.
  `endTime 0.06`, fixed dt = 1e-5 → 6000 steps total; window =
  `nStepsFrac = 0.6667` selects the **last 4000 samples** (t = 0.02001…0.06 s).
- **Tier 2 (confirmation): t0 = 0.100 s, t1 = 0.300 s** (case endTime),
  T_win = 200 ms ≈ 13.8 transits, `endTime 0.3`,
  `nStepsFrac = 2/3` → last **20000 samples** (t = 0.10001…0.3 s), 30000
  steps total. Run only after Tier 1 passes.

Stationarity gate before using a window: mean of the first and second half
of the baseline window must differ by < 1 σ_J; otherwise move t0 later.

## 6. Finite-difference protocol

- **Central differences:** dJ/dα ≈ [J(1+h) − J(1−h)] / (2h), with
  **h ∈ {0.005, 0.01, 0.02}** (0.5 %, 1 %, 2 % of the jet velocity).
- **Curvature diagnostic:** C = [J(1+h) + J(1−h) − 2J0] / h² — reported for
  each h; a large C relative to the first-order term means h is too large
  (reduce h) or the response is genuinely nonlinear (report both).
- **SNR rule for the usable h:** accept an h only if
  **|J(1+h) − J(1−h)| ≥ 5·√2·σ_J** (§7); take the smallest h satisfying it,
  and require the estimates for the two largest valid h's to agree within
  25 %.
- **Pilot step-size probe:** before the full matrix, run J0 and one
  ±1 % pair; if |dJ/dα| is below the noise floor at h = 1 %, escalate to
  h = 2 % (and, if still marginal, widen the window rather than the step).

## 7. Repeatability checks and sampling uncertainty

**Repeatability (all must pass before interpreting any FD number):**

- **R1 determinism:** run the baseline twice from identical inputs →
  per-step `meanTOutlet:` series and final J must be **bit-identical**
  (max abs diff = 0). Non-zero difference means the platform/binary is not
  reproducible (parallel reduction order, BC recompile variance) → stop and
  investigate; FD is meaningless until it passes.
- **R2 dt sequence:** the `Time =` series must be identical across all runs
  (fixed dt makes this trivial; assert it anyway — it also catches a
  silently re-enabled `adjustTimeStep`).
- **R3 objective plumbing:** recompute J offline from the per-step log values
  and require agreement with `evalFunctions` to ≤1e-9 relative (as in D2,
  7.6e-12); with fixed dt the unweighted and dt-weighted means must be equal.
- **R4 window membership:** number of samples in [t0, t1] must equal
  T_win/dt exactly.

**Sampling uncertainty.** For a window of length T_win with per-sample σ_inst
and integral time scale τ_int (measured *in the baseline window*, not from
the early D2 data):

    σ_J = σ_inst · sqrt(2·τ_int / T_win)          (effective independent
                                                   samples N_eff = T_win/(2τ_int))
    FD difference noise      = √2·σ_J
    95 % CI on dJ/dα         = ±1.96·√2·σ_J / (2h)

Measured anchors (early window, will be superseded by the baseline
measurement): σ_inst = 1.28e-2 K, τ_int ≈ 0.87 ms. Illustration with
τ_int = 5 ms, σ_inst = 0.05 K, T_win = 40 ms → σ_J = 0.025 K, so at
h = 0.01 a detectable sensitivity must exceed
1.96·√2·0.025/0.01 ≈ **7 K per unit α**; a response of order 1 K per unit α
would need T_win ≥ 400 ms (Tier 2) or h = 2 %. This is why σ_J is measured
first and h is chosen by the SNR rule rather than fixed up front.

**Note on chaos:** σ_J falls only as 1/√N_eff; extending T_win past a few
dozen τ_int buys little if the finite-time Lyapunov growth keeps the two
±h trajectories decorrelated — Tier 2 exists to check whether the estimate
actually stabilizes with window length (that *is* the identifiability test).

## 8. Acceptance criteria

| id | criterion | threshold |
|---|---|---|
| G1 | determinism (R1) | bit-identical baseline repeat |
| G2 | identical dt sequence (R2) | `Time =` series equal across runs |
| G3 | plumbing (R3/R4) | J match ≤1e-9 rel; sample count exact |
| G4 | stationarity of [t0, t1] | half-window means differ < 1 σ_J |
| G5 | signal-to-noise | \|J+ − J−\| ≥ 5·√2·σ_J for the reported h |
| G6 | step-size consistency | central-difference estimates for the two largest valid h agree within 25 %, same sign |
| G7 | primal health | `primalFail = 0`, no NaN/Inf in the per-step series, Co ≤ 0.4 in the pilot |

Report: dJ/dα with the 95 % CI from §7, the curvature C, the window, and
whether the effect is distinguishable from zero (CI excludes 0). A result
failing G5 is reported as **not identifiable at this window length**, not as
a sensitivity value.

## 9. Identifiability assessment (what this study can and cannot conclude)

- *Structural:* one scalar control → one scalar objective; J(α) is
  identifiable in principle whenever ∂J/∂α ≠ 0 (no parameter sloppiness,
  no null space).
- *Practical:* limited by (i) chaotic sampling noise (σ_J, §7) — N_eff
  independent samples per window, (ii) the curvature/h trade-off, and
  (iii) dt-path contamination if fixed dt is violated. The planned checks
  (G5/G6, Tier 1 → Tier 2 window extension) separate these: G6-failing =
  truncation error, Tier-1/2 disagreement = sampling error.
- *Known-weak control:* ambient H2O is structurally inert in the current
  configuration (§2) — defer, do not read a null result as noise.

## 10. Run matrix and cost

Measured rate: 0.86 s/step serial (D2 run: 229.85 s / 266 steps).

| tier | runs | steps/run | wall time/run | total |
|---|---|---|---|---|
| pilot | 1 baseline at fixed dt to t1 (Co check) | 6000 | ≈ 86 min | 1.4 h |
| Tier 1 screening | J0 ×2 (R1), ±h for h ∈ {0.01, 0.02} → 6 | 6000 | ≈ 86 min | ≈ 8.6 h |
| Tier 1 full | + ±h at h = 0.005 → 8 runs | 6000 | ≈ 86 min | ≈ 11.5 h |
| Tier 2 confirmation | J0, ±h at best h → 3 | 30000 | ≈ 7.2 h | ≈ 21.5 h |

(A larger Co-safe fixed dt found in the pilot divides these numbers
proportionally. Runs can be overlapped across cores if more than one
process is allowed; estimates assume the current single-serial run setup.)

## 11. Deliverables / next steps

1. Pilot: fixed-dt stability + stationarity of the Tier-1 window; measure
   σ_inst, τ_int, σ_J in that window.
2. R1 determinism gate; then the Tier-1 FD matrix; report dJ/dα with CI
   under G1–G7.
3. Tier-2 confirmation if Tier 1 passes G5–G6.
4. Feed the *same* short-window objective + control to Phase D4 as the FD
   reference for ADR validation (D4 is the audit; its rebuild/FD gate stays
   closed until explicitly requested).
5. Deferred: ambient-H2O control (post-Phase E), dt-weighted timeOp if
   exact `fieldAverage` parity is ever required (D1 §8).
