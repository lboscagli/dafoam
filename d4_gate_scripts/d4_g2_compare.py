#!/usr/bin/env python3
"""Phase D4 gate G2 verdict: compare the two residual dumps (<=1e-10).

    d4_g2_compare.py g2_residuals_orig.npy g2_residuals_adr.npy
"""
import sys

import numpy as np

TOL = 1e-10


def main():
    orig = np.load(sys.argv[1]).astype(np.float64)
    adr = np.load(sys.argv[2]).astype(np.float64)

    print("G2: orig n=%d  adr n=%d" % (orig.size, adr.size), flush=True)
    if orig.shape != adr.shape:
        print("G2: FAIL (residual vector sizes differ)", flush=True)
        return 2

    scale = float(np.max(np.abs(orig)))
    if scale == 0.0:
        scale = float(np.max(np.abs(adr)))
    absdiff = np.abs(orig - adr)
    max_abs = float(absdiff.max())
    max_rel = float(absdiff.max() / scale) if scale else 0.0

    nz = np.abs(orig) > 0.0
    per_elem = np.zeros_like(absdiff)
    per_elem[nz] = absdiff[nz] / np.abs(orig[nz])
    max_per_elem = float(per_elem.max()) if nz.any() else 0.0

    print("G2: scale=%.17g" % scale, flush=True)
    print("G2: max_abs_diff=%.17g" % max_abs, flush=True)
    print("G2: max_rel_diff=%.17g  (max elementwise rel=%.17g)"
          % (max_rel, max_per_elem), flush=True)
    print("G2: elements differing at all: %d / %d"
          % (int((absdiff > 0.0).sum()), absdiff.size), flush=True)
    print("G2: elements beyond 1e-10 (scaled): %d"
          % int((per_elem > TOL).sum()), flush=True)

    verdict = max_per_elem <= TOL
    print("G2: %s (max elementwise rel %.3g <= %.0e)"
          % ("PASS" if verdict else "FAIL", max_per_elem, TOL), flush=True)
    return 0 if verdict else 3


if __name__ == "__main__":
    sys.exit(main())
