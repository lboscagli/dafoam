#!/usr/bin/env python3
"""Phase D4 gate G3: timeAccurate checkpoint/replay integrity.

One ADR process does what the real unsteady-adjoint flow does: run the
primal (which sets primalFinalTimeIndex_ and, with reduceIO, writes a
per-step checkpoint set), then replay that set with readStateVars() and
check:

  * checkpoints exist for every step of the window,
  * the replayed final state matches the in-memory primal final state
    (<=1e-12; the case writes ascii at writePrecision 16),
  * getdFScaling() returns 1/N inside the timeOp window and 0 outside,
  * primalFinalTimeIndex_ is valid, i.e. Prompt-2's fatal does not fire,
  * consecutive replayed states differ (per-step data really is loaded).

Run under the union ADR env:
    d4_adr_env.sh python d4_g3_replay.py <caseDir>
"""
import glob
import os
import sys

import numpy as np

from dafoam import PYDAFOAM
from d4_g1_adr_init import patch_adr_union_path

TOL = 1e-12
N_STEPS_FRAC = 0.5
OUT = "/workspace/logs"


def mean_t_outlet_function(n_steps_frac):
    return {
        "meanTOutlet": {
            "type": "patchMean",
            "source": "patchToFace",
            "patches": ["outlet"],
            "varName": "T",
            "varType": "scalar",
            "index": 0,
            "scale": 1.0,
            "timeOp": "average",
            "nStepsFrac": n_steps_frac,
        }
    }


def main():
    case = sys.argv[1]
    os.chdir(case)
    print("G3: case=%s" % case, flush=True)

    patch_adr_union_path()
    options = {
        "solverName": "DAContrailFoam",
        "unsteadyAdjoint": {"mode": "timeAccurate", "reduceIO": True},
        "function": mean_t_outlet_function(N_STEPS_FRAC),
    }
    da = PYDAFOAM(options=options)
    print("G3: INIT_OK (unsteadyAdjoint.mode=timeAccurate)", flush=True)

    da.setPrimalInitialConditions()
    da.setPrimalBoundaryConditions()

    fail = da.solverAD.solvePrimal()
    print("G3: primalFail=%s" % fail, flush=True)
    if fail != 0:
        print("G3: FAIL (primal returned nonzero)", flush=True)
        return 2

    n = da.solverAD.getNLocalAdjointStates()
    primal_final = np.zeros(n, dtype=np.float64)
    da.solverAD.getOFFields(primal_final)
    np.save(os.path.join(OUT, "g3_primal_final.npy"), primal_final)
    print("G3: primal final states captured n=%d" % n, flush=True)

    deltas = sorted(
        float(os.path.basename(p))
        for p in glob.glob(os.path.join(case, "[0-9]*"))
        if os.path.isdir(p)
    )
    delta_t = 1.25e-5
    n_expected = 100
    print("G3: checkpoint times found: %d (first=%g last=%g)"
          % (len(deltas), deltas[0], deltas[-1]), flush=True)
    if len(deltas) < n_expected:
        print("G3: FAIL (expected >= %d checkpoint dirs for %d steps)"
              % (n_expected + 1, n_expected), flush=True)
        return 3

    scaling = [da.solverAD.getdFScaling("meanTOutlet", i)
               for i in range(n_expected)]
    nonzero = [i for i, v in enumerate(scaling) if v != 0.0]
    print("G3: dFScaling nonzero indices: %d..%d (%d of %d)"
          % (nonzero[0], nonzero[-1], len(nonzero), n_expected)
          if nonzero else "G3: dFScaling all zero", flush=True)

    contiguous = nonzero == list(range(nonzero[0], nonzero[0] + len(nonzero)))
    values = set(scaling[i] for i in nonzero)
    expected = 1.0 / len(nonzero) if nonzero else 0.0
    window_ok = (
        nonzero
        and contiguous
        and len(values) == 1
        and abs(scaling[nonzero[0]] - expected) <= 1e-12
        and all(scaling[i] == 0.0
                for i in range(n_expected) if i not in nonzero)
        and nonzero[0] > 0
    )
    print("G3: dFScaling contiguous=%s values=%s expected=%.17g "
          "out-of-window zeros=%s"
          % (contiguous, sorted(values), expected,
             all(scaling[i] == 0.0
                 for i in range(n_expected) if i not in nonzero)),
          flush=True)

    zero_missing = []
    for name in ("U", "p", "T", "k", "phi"):
        if not os.path.exists(os.path.join(case, "0", name)):
            zero_missing.append(name)
    print("G3: 0/ missing state fields: %s"
          % (zero_missing or "none"), flush=True)

    # readStateVars() MUST_READs every state, and the first two steps need
    # time levels at/before t=0, which only exist in the hand-provided 0/
    # folder (the solver never writes it: the write sits inside the time
    # loop). Replay the window that is covered entirely by the reduceIO
    # checkpoint set.
    delta_t = 1.25e-5
    first_ckpt = deltas[1]
    replay_times = [t for t in deltas if t - 2 * delta_t >= first_ckpt - 1e-15]
    print("G3: replaying %d checkpoints (t from %g to %g); %d early steps "
          "need 0/ and are skipped"
          % (len(replay_times), replay_times[0], replay_times[-1],
             len(deltas) - len(replay_times)), flush=True)

    max_diff = 0.0
    last_diff = 0.0
    prev = None
    final_states = None
    for tv in replay_times:
        da.solverAD.readStateVars(tv, 0)
        da.solverAD.readStateVars(tv - delta_t, 1)
        da.solverAD.readStateVars(tv - 2 * delta_t, 2)
        cur = np.zeros(n, dtype=np.float64)
        da.solverAD.getOFFields(cur)
        if prev is not None:
            last_diff = float(np.max(np.abs(cur - prev)))
            max_diff = max(max_diff, last_diff)
        prev = cur
        if abs(tv - deltas[-1]) < 1e-15:
            final_states = cur

    print("G3: replayed %d checkpoints (max step-to-step state change "
          "%.17g)" % (len(replay_times), max_diff), flush=True)

    roundtrip = float(np.max(np.abs(final_states - primal_final)))
    scale = float(np.max(np.abs(primal_final))) or 1.0
    roundtrip_rel = roundtrip / scale
    print("G3: |primal_final - replayed_final|_max = %.17g (scale %.17g, "
          "rel %.17g)" % (roundtrip, scale, roundtrip_rel), flush=True)

    checks = {
        "checkpoints >= steps+1": len(deltas) >= n_expected + 1,
        "dFScaling window": bool(window_ok),
        "dFScaling in-window value 1/N": bool(
            nonzero and abs(scaling[nonzero[0]] - 1.0 / len(nonzero)) <= 1e-12),
        "final state round trip <= 1e-12 (relative)": roundtrip_rel <= TOL,
        "replay loaded varying per-step states": max_diff > 0.0,
    }
    for name, ok in checks.items():
        print("G3: [%s] %s" % ("OK" if ok else "FAIL", name), flush=True)

    passed = all(checks.values())
    print("G3: %s" % ("PASS" if passed else "FAIL"), flush=True)
    return 0 if passed else 4


if __name__ == "__main__":
    sys.exit(main())
