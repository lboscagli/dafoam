#!/usr/bin/env bash
# Build the two OpenFOAM functionObject libs that the D4 case's controlDict
# loads - "coded" (libutilityFunctionObjects) and "fieldAverage"
# (libfieldFunctionObjects) - in the active AD mode (WM_AD_MODE=ADR/ADF).
#
# Why this is needed (Phase D4 / gate G1):
#   The AD OpenFOAM tree ships NO functionObject libs at all (96 libs in
#   .../linux64GccDPInt32OptADR/lib, zero from src/functionObjects). The ADR
#   process therefore dlopened the *unsuffixed* (double-built) libs from the
#   v2506 tree via the case's controlDict `libs` entries -> the double
#   codedFunctionObject then loaded an AD/codi-generated dynamicCode lib ->
#   ABI mismatch (instant/scalar layout) -> SEGV during solvePrimal.
#
# Pre-reqs done in the same session:
#   - src/functionObjects/{utilities,field}/Make/options patched so every -l
#     in LIB_LIBS carries $(WM_AD_MODE) (option B: source patch, same as
#     src/finiteVolume/Make/options). Backup: Make/options.orig-d4
#   - codedFunctionObject.C patched in BOTH the AD and the v2506 trees
#     (backup: codedFunctionObject.C.orig-d4)
#
# Naming: wmake emits LIB = $(FOAM_LIBBIN)/lib<name> (unsuffixed), and the
# case's controlDict requests the unsuffixed names. In the ADR union env the
# AD platform lib dir precedes the v2506 dir in LD_LIBRARY_PATH, so the AD
# copy wins; in the original env only the v2506 dir is on the path, so the
# double copy is used. -> the SAME case works in both modes, no dict edit.
#
# Usage (inside dafoam-dev):
#   d4_adr_env.sh bash /workspace/scripts/d4_build_functionobjects.sh
set -o pipefail

readonly AD_ROOT=/home/dafoamuser/dafoam/OpenFOAM/OpenFOAM-AD
readonly SHIM_DIR=/tmp/wmake-shim
readonly NCOMP=2

[ -n "$WM_AD_MODE" ] || { echo "FATAL: no AD env (run via d4_adr_env.sh)" >&2; exit 3; }

export WM_NCOMPPROCS="$NCOMP"
# OF bashrc prepends its own bin dirs -> force the bare -j -> -j2 shim first
filtered=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vxF "$SHIM_DIR" | paste -sd:)
export PATH="$SHIM_DIR:$filtered"
hash -r
echo "wmake -> $(command -v wmake)  WM_AD_MODE=$WM_AD_MODE  WM_OPTIONS=$WM_OPTIONS"

cd "$AD_ROOT" || exit 1

for d in src/functionObjects/utilities src/functionObjects/field; do
    echo "=== wmakeLnInclude $d  $(date -u -Iseconds)"
    wmakeLnInclude "$d" || { echo "FAILED lnInclude $d" >&2; exit 1; }
done

for d in src/functionObjects/utilities src/functionObjects/field; do
    echo "=== wmake $d  $(date -u -Iseconds)"
    ( cd "$d" && wmake -j"$NCOMP" libso ) || { echo "FAILED build $d" >&2; exit 1; }
    echo "=== done $d  $(date -u -Iseconds)"
done

LIBDIR=$AD_ROOT/platforms/$WM_OPTIONS/lib
echo "=== outputs in $LIBDIR:"
ls -la "$LIBDIR"/libutilityFunctionObjects.so "$LIBDIR"/libfieldFunctionObjects.so || exit 1
echo "=== codedFunctionObject present in libutilityFunctionObjects.so:"
strings "$LIBDIR"/libutilityFunctionObjects.so | grep -c codedFunctionObject
echo "=== NEEDED of both libs (all must be *$WM_AD_MODE):"
for l in libutilityFunctionObjects.so libfieldFunctionObjects.so; do
    echo "-- $l"
    readelf -d "$LIBDIR/$l" | grep NEEDED
    bad=$(readelf -d "$LIBDIR/$l" | grep NEEDED | grep -vc "${WM_AD_MODE}\.so")
    echo "   unsuffixed-NEEDED count = $bad (want 0)"
done
echo "=== build completed: $(date -u -Iseconds)"
