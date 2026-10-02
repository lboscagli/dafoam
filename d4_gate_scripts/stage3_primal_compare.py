#!/usr/bin/env python
"""Stage 3 verification: DAContrailFoam primal vs standalone contrailFoam.

Modes
-----
--mode da  --case DIR --out FILE.npz
    Build DAContrailFoam, run solvePrimal on DIR, dump the final-time fields.

--mode of  --case DIR --out FILE.npz
    Build DAContrailFoam (no solve) on a finished contrailFoam case so the
    fields of its latest time directory are read through the same public API.

--mode cmp --a FILE.npz --b FILE.npz
    Compare two dumps cell by cell (b is the reference).

--func (mode da only)
    Inject the meanTOutlet objective (area-weighted patch mean of T on the
    outlet patch, timeOp average over the full run window) into daOptions and
    print the evalFunctions result after solvePrimal. With nStepsFrac=1.0 the
    time-op window covers [0, endTime], i.e. every recorded step.

Bypasses ADR/ADF initialization exactly like smoke_test_dacontrail.py: the AD
libraries are stale and COMPILE_DAFOAM_NOAD skips rebuilding them.
"""

import argparse
import os
import sys

import numpy as np


SCALARS = ["p", "T", "O2", "CO2", "H2O", "k", "nut", "rho"]
VECTORS = ["U"]


def _initSolver_original_only(self):
    if self.solverInitialized == 1:
        raise RuntimeError("_initSolver called twice")

    solverName = self.getOption("solverName")
    solverArg = solverName + " -python " + self.parallelFlag

    from dafoam.libs.pyDASolvers import pyDASolvers

    self.solver = pyDASolvers(solverArg.encode(), self.options)
    self.solverAD = self.solver
    self.solver.initSolver()
    self.solverInitialized = 1


def _mean_t_outlet_function(n_steps_frac=1.0):
    """meanTOutlet: area-weighted mean of T over the outlet patch, averaged in
    time over a trailing nStepsFrac window (1.0 = the whole run window
    [0, endTime]; (t1-t0)/t1 selects [t0, t1] with endTime = t1)."""
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


def _patch_and_build(case, with_func=False, n_steps_frac=1.0):
    os.chdir(case)

    import dafoam.pyDAFoam as pyDAFoamModule
    from dafoam import PYDAFOAM

    pyDAFoamModule.PYDAFOAM._initSolver = _initSolver_original_only

    daOptions = {"solverName": "DAContrailFoam"}
    if with_func:
        daOptions["function"] = _mean_t_outlet_function(n_steps_frac)
    return PYDAFOAM(options=daOptions)


def dump(dafoam, path):
    n_cells = dafoam.solver.getNLocalCells()
    out = {"nCells": np.array([n_cells], dtype=np.int64)}
    for name in SCALARS:
        arr = np.zeros(n_cells, dtype="d")
        dafoam.solver.getOFField(name, "scalar", arr)
        out[name] = arr
    for name in VECTORS:
        arr = np.zeros(3 * n_cells, dtype="d")
        dafoam.solver.getOFField(name, "vector", arr)
        out[name] = arr
    np.savez(path, **out)
    print("STAGE3: dumped %d cells to %s" % (n_cells, path), flush=True)


def compare(a_path, b_path, rtol):
    a = np.load(a_path)
    b = np.load(b_path)

    if int(a["nCells"][0]) != int(b["nCells"][0]):
        print("STAGE3: FAIL cell count %d vs %d"
              % (int(a["nCells"][0]), int(b["nCells"][0])))
        return 1

    ok = True
    print("STAGE3: %-10s %14s %14s %14s %10s"
          % ("field", "maxAbsDiff", "maxAbsRef", "maxRelDiff", "verdict"))
    for name in SCALARS + VECTORS:
        da = a[name]
        db = b[name]
        max_abs = float(np.max(np.abs(da - db)))
        max_ref = float(np.max(np.abs(db)))
        max_rel = max_abs / max(max_ref, 1e-300)
        verdict = "OK" if max_rel <= rtol else "DIFF"
        if max_rel > rtol:
            ok = False
        print("STAGE3: %-10s %14.6e %14.6e %14.6e %10s"
              % (name, max_abs, max_ref, max_rel, verdict))

    # k is a small part of the cell energy; also judge its error against the
    # total kinetic energy TKE = k + 0.5*|U|^2 of the reference (b), both
    # cell-wise (max) and energy-weighted over the whole domain (sum).
    ub = b[VECTORS[0]].reshape(-1, 3)
    tke = b["k"] + 0.5 * np.einsum("ij,ij->i", ub, ub)
    dk = np.abs(a["k"] - b["k"])
    for tag, ratio in (
        ("k/TKE", dk / np.maximum(tke, 1e-300)),
        ("k/TKEsum", None),
    ):
        if ratio is None:
            rel = float(dk.sum() / max(tke.sum(), 1e-300))
            max_abs = float(dk.sum())
            max_ref = float(tke.sum())
        else:
            rel = float(np.max(ratio))
            i = int(np.argmax(ratio))
            max_abs = float(dk[i])
            max_ref = float(tke[i])
        verdict = "OK" if rel <= rtol else "DIFF"
        if rel > rtol:
            ok = False
        print("STAGE3: %-10s %14.6e %14.6e %14.6e %10s"
              % (tag, max_abs, max_ref, rel, verdict))

    print("STAGE3 RESULT: %s" % ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=["da", "of", "cmp"])
    parser.add_argument("--case")
    parser.add_argument("--out")
    parser.add_argument("--a")
    parser.add_argument("--b")
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument(
        "--func",
        action="store_true",
        help="inject the meanTOutlet patchMean(T,outlet) timeOp average "
        "function dict and print evalFunctions after the primal run",
    )
    parser.add_argument(
        "--nStepsFrac",
        type=float,
        default=1.0,
        help="trailing window fraction for --func (1.0 = [0, endTime]; "
        "(t1-t0)/t1 selects [t0, t1] when endTime = t1)",
    )
    args = parser.parse_args()

    if args.mode == "cmp":
        raise SystemExit(compare(args.a, args.b, args.rtol))

    if args.func and args.mode != "da":
        raise SystemExit("--func requires --mode da (a primal run to average)")

    if args.nStepsFrac != 1.0 and not args.func:
        raise SystemExit("--nStepsFrac requires --func")

    dafoam = _patch_and_build(args.case, with_func=args.func,
                              n_steps_frac=args.nStepsFrac)

    if args.mode == "da":
        print("STAGE3: running DAContrailFoam.solvePrimal", flush=True)
        dafoam()
        print("STAGE3: solvePrimal done (primalFail=%d)" % dafoam.primalFail,
              flush=True)
        if dafoam.primalFail != 0:
            raise SystemExit(2)

        if args.func:
            print("STAGE3: func window nStepsFrac=%.6f" % args.nStepsFrac,
                  flush=True)
            funcs = {}
            dafoam.evalFunctions(funcs)
            for funcName in sorted(funcs):
                print("STAGE3: func %s = %.12e"
                      % (funcName, float(funcs[funcName])), flush=True)

    dump(dafoam, args.out)
    raise SystemExit(0)


if __name__ == "__main__":
    main()
