#!/usr/bin/env bash
# Original (double) runtime env for D4 cross-mode gate runs.
# Source chain stops at loadDAFoam.sh: no OpenFOAM-AD bashrc, so WM_AD_MODE /
# WM_OPTIONS stay on the v2506 platform and only the unsuffixed libs resolve.
# Counterpart of d4_adr_env.sh (q.v. for why "$@" is cleared first).
CMD=("$@")
set --
set +e
source /home/dafoamuser/dafoam/loadDAFoam.sh
set -e
[ -z "${WM_AD_MODE:-}" ] || { echo "FATAL: WM_AD_MODE=$WM_AD_MODE set in original env" >&2; exit 3; }
exec "${CMD[@]}"
