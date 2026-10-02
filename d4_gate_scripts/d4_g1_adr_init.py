#!/usr/bin/env python3
"""Phase D4 gate G1: ADR-mode initSolver + short ADR primal on a case.

Answers audit section 4 question 1: does initSolver() + solvePrimal() run in
ADR mode, including the runtime-compiled coded BCs (dynamicCode)?

Run under the union ADR env:
    d4_adr_env.sh python d4_g1_adr_init.py <caseDir>
Expects a short case (e.g. endTime 0.00125 = 100 steps at deltaT 1.25e-5)
without a dynamicCode/ dir so the BC compile happens fresh under ADR.

G1 runs the union path (PYDAFOAM builds both the original and the ADR
solver bindings) but the ADR init sequence differs from the plain double
path in three ways, each recorded in the runbook G1 fix table: the original
solver is constructed but not initialized (its initSolver needs the
non-AD scalar ABI), pyofm's mesh read has to happen before the ADR
initSolver, and the setPrimal* helpers are ADR-only.
"""
import os
import sys
import time

from dafoam import PYDAFOAM


def patch_adr_union_path():
    orig_read_mesh_info = PYDAFOAM._readMeshInfo

    def _readMeshInfo(self):
        if getattr(self, "_meshInfoDone", False):
            print("G1: _readMeshInfo cached (already ran before ADR initSolver)",
                  flush=True)
            return
        orig_read_mesh_info(self)
        self._meshInfoDone = True

    def _initSolver(self):
        solverName = self.getOption("solverName")
        solverArg = solverName + " -python " + self.parallelFlag
        from dafoam.libs.pyDASolvers import pyDASolvers
        from dafoam.libs.ADR.pyDASolvers import pyDASolvers as pyDASolversAD
        self.solver = pyDASolvers(solverArg.encode(), self.options)
        print("G1: original solver bindings constructed (initSolver skipped)",
              flush=True)
        self.solverAD = pyDASolversAD(solverArg.encode(), self.options)
        print("G1: ADR solver bindings constructed", flush=True)
        self._readMeshInfo()
        print("G1: pyofm mesh read done BEFORE ADR initSolver", flush=True)
        self.solverAD.initSolver()
        print("G1: ADR initSolver done", flush=True)
        self.solverInitialized = 1

    def _setPrimalInitialConditions(self, printInfo=1, printInfoAD=0):
        self.solverAD.setPrimalInitialConditions(printInfoAD)
        self.solverAD.getInitStateVals(printInfoAD)

    def _setPrimalBoundaryConditions(self, printInfo=1, printInfoAD=0):
        self.solverAD.setPrimalBoundaryConditions(printInfoAD)

    PYDAFOAM._readMeshInfo = _readMeshInfo
    PYDAFOAM._initSolver = _initSolver
    PYDAFOAM.setPrimalInitialConditions = _setPrimalInitialConditions
    PYDAFOAM.setPrimalBoundaryConditions = _setPrimalBoundaryConditions


def main():
    case = sys.argv[1]
    os.chdir(case)
    print("G1: case=%s cwd=%s" % (case, os.getcwd()), flush=True)

    patch_adr_union_path()

    t0 = time.time()
    da = PYDAFOAM(options={"solverName": "DAContrailFoam"})
    print("G1: INIT_OK elapsed=%.1fs (ADR solver initialized, "
          "original solver bindings loaded)" % (time.time() - t0), flush=True)

    t0 = time.time()
    fail = da.solverAD.solvePrimal()
    elapsed = time.time() - t0
    print("G1: ADR_PRIMAL primalFail=%s elapsed=%.1fs" % (fail, elapsed),
          flush=True)
    if fail != 0:
        print("G1: FAIL (ADR primal returned nonzero)", flush=True)
        return 2

    print("G1: PASS (ADR init + primal + coded-BC compile OK)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
