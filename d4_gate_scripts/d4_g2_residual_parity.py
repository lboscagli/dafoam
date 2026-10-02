#!/usr/bin/env python3
"""Phase D4 gate G2: residual parity, original mode vs ADR mode.

Dumps the full-precision residual vector for one mode on the reference case
so the two runs can be compared offline (audit section 9: agreement <=1e-10
on the same case/states). Uses getResiduals() rather than
calcPrimalResidualStatistics("print") because OF's Info stream only carries
6 significant digits.

Run once per mode, two processes / two envs, same case directory:
    d4_orig_env.sh python d4_g2_residual_parity.py orig /tmp/stage3/d4_g1
    d4_adr_env.sh  python d4_g2_residual_parity.py adr  /tmp/stage3/d4_g1
then:
    d4_g2_compare.py g2_residuals_orig.npy g2_residuals_adr.npy
"""
import os
import sys

import numpy as np

from dafoam import PYDAFOAM
from d4_g1_adr_init import patch_adr_union_path


def patch_original_path():
    def _initSolver(self):
        solverName = self.getOption("solverName")
        solverArg = solverName + " -python " + self.parallelFlag
        from dafoam.libs.pyDASolvers import pyDASolvers
        self.solver = pyDASolvers(solverArg.encode(), self.options)
        self.solverAD = self.solver
        self.solver.initSolver()
        self.solverInitialized = 1

    def _setPrimalInitialConditions(self, printInfo=1, printInfoAD=0):
        self.solver.setPrimalInitialConditions(printInfo)

    def _setPrimalBoundaryConditions(self, printInfo=1, printInfoAD=0):
        self.solver.setPrimalBoundaryConditions(printInfo)

    PYDAFOAM._initSolver = _initSolver
    PYDAFOAM.setPrimalInitialConditions = _setPrimalInitialConditions
    PYDAFOAM.setPrimalBoundaryConditions = _setPrimalBoundaryConditions


def main():
    mode = sys.argv[1]
    case = sys.argv[2]
    os.chdir(case)
    print("G2: mode=%s case=%s" % (mode, case), flush=True)

    if mode == "adr":
        patch_adr_union_path()
    elif mode == "orig":
        patch_original_path()
    else:
        print("G2: FAIL (unknown mode %r)" % mode, flush=True)
        return 2

    da = PYDAFOAM(options={"solverName": "DAContrailFoam"})
    print("G2: mode=%s INIT_OK" % mode, flush=True)

    da.setPrimalInitialConditions()
    da.setPrimalBoundaryConditions()
    print("G2: primal initial + boundary conditions set", flush=True)

    solver = da.solverAD if mode == "adr" else da.solver
    n = solver.getNLocalAdjointStates()
    residuals = np.zeros(n, dtype=np.float64)
    solver.getResiduals(residuals)

    out = os.path.join("/workspace/logs", "g2_residuals_%s.npy" % mode)
    np.save(out, residuals)
    print("G2: n=%d out=%s" % (n, out), flush=True)
    print("G2: residual_norm2=%.17g absmax=%.17g"
          % (float(np.linalg.norm(residuals)),
             float(np.max(np.abs(residuals)))), flush=True)
    print("G2: DUMP_OK mode=%s" % mode, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
