#!/usr/bin/env bash
# Rebuild OpenFOAM-AD libOpenFOAM in the active AD mode (WM_AD_MODE set by
# d4_adr_env.sh). Needed after the scalarField.{H,C} sumProd<scalar>
# specialisation: the generic sumProd uses OF's "dotdot" (&&) operator, which
# is correct for vector/tensor but yields a logical-AND count for scalar-like
# types; double/float mask it via their * specialisations, codi ActiveType
# (scalar under ADR/ADF) did not - so gSumProd was wrong for every BiCG/CG
# solver. Objects are not cached in this image, so this is a full build.
# Run: D4_AD_MODE=ADR d4_adr_env.sh bash /workspace/scripts/d4_rebuild_libOpenFOAM.sh
set -euo pipefail
cd "$FOAM_SRC/OpenFOAM"
echo "WM_OPTIONS=$WM_OPTIONS  WM_AD_MODE=$WM_AD_MODE  FOAM_LIBBIN=$FOAM_LIBBIN"
date
wmake -j2 libso
date
echo "BUILD_DONE"
