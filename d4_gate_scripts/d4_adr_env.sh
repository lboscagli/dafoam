#!/usr/bin/env bash
# Union ADR runtime env for DAContrailFoam Phase D4 gate runs.
# - loadDAFoam.sh: original env, so the original-mode solver that pyDAFoam's
#   _initSolver constructs alongside the ADR one resolves its non-suffixed libs
#   (verified: original pyDASolvers .so has 37 not-found deps under pure ADR env)
# - OpenFOAM-AD bashrc with WM_AD_MODE=ADR: AD headers/wmake for the
#   runtime-compiled coded BCs (audit section 4) and the *ADR-suffixed libs
#   (verified: ADR pyDASolvers .so has 0 not-found deps under this env)
# - original platform lib dirs appended AFTER the AD ones: SONAMEs are disjoint
#   (libOpenFOAM.so vs libOpenFOAMADR.so) so both can coexist in one process;
#   libPstream/libptscotchDecomp live in lib/sys-openmpi
# Usage: d4_adr_env.sh <cmd> [args...]
# D4_AD_MODE=ADR|ADF selects the AD mode (default ADR).
# Save the command first and clear "$": some env script in the source chain
# iterates "$@" and dot-sources/evals file-like args (observed: it ran
# `d4_g1_adr_init.py` as a bash script because it was arg #3 of this wrapper).
CMD=("$@")
set --
set -e
D4_AD_MODE=${D4_AD_MODE:-ADR}
case "$D4_AD_MODE" in ADR|ADF) ;; *) echo "FATAL: bad D4_AD_MODE=$D4_AD_MODE" >&2; exit 4;; esac
DEP_ROOT=/home/dafoamuser/dafoam
AD_ROOT=$DEP_ROOT/OpenFOAM/OpenFOAM-AD
OF_ROOT=$DEP_ROOT/OpenFOAM/OpenFOAM-v2506
OF_LIB=$OF_ROOT/platforms/linux64GccDPInt32Opt/lib

# OpenFOAM's config.sh/setup errors under set -e (pop_var_context) -> errexit
# off while sourcing env scripts, on otherwise.
set +e
source "$DEP_ROOT/loadDAFoam.sh"
set -e
sed -i "s/export WM_AD_MODE=.*/export WM_AD_MODE=$D4_AD_MODE/" "$AD_ROOT/etc/bashrc"
set +e
source "$AD_ROOT/etc/bashrc"
set -e
[ "$WM_AD_MODE" = "$D4_AD_MODE" ] || { echo "FATAL: $D4_AD_MODE env not active (got WM_AD_MODE=$WM_AD_MODE)" >&2; exit 3; }
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$OF_LIB:$OF_LIB/sys-openmpi"
exec "${CMD[@]}"
