#!/usr/bin/env bash
# OOM-safe ADR + ADF rebuild for DAContrailFoam (Phase D4 gate G1).
# Run inside dafoam-dev:
#   docker exec -d dafoam-dev bash --noprofile --norc /workspace/scripts/build_ad_adf.sh
# Differences vs src/dafoam/Allmake:
#  - strict set -e with subshell cd's (Allmake's &&-chains mask wmake failures
#    and desync the cwd, so a failed pass still reports SUCCEEDED)
#  - wmake job shim: OpenFOAM's config.sh/unset clears WM_NCOMPPROCS and bare
#    `wmake -j` then resets to nproc=10 (env ignored) -> OOM-kills cc1plus in
#    the 7.8 GiB Docker VM on codi-heavy AD templates. The shim rewrites bare
#    `-j` to `-j2` for every wmake invocation.
set -o pipefail

readonly DEP_ROOT=/home/dafoamuser/dafoam
readonly SOURCE_DIR=/workspace/src/dafoam
readonly LOG_DIR=/workspace/logs
readonly STAMP="$(date -u +%Y%m%d-%H%M%S)"
readonly LOG_FILE="$LOG_DIR/dafoam-Allmake-$STAMP.log"
readonly STATUS_FILE="$LOG_DIR/dafoam-Allmake-status.txt"
readonly PID_FILE="$LOG_DIR/dafoam-Allmake.pid"
readonly LATEST_FILE="$LOG_DIR/dafoam-Allmake-latest.log"
readonly NCOMP=2

mkdir -p "$LOG_DIR"
printf '%s\n' "$LOG_FILE" > "$LATEST_FILE"
printf '%s\n' "$$" > "$PID_FILE"
printf 'STARTED time=%s\nlog=%s\n' "$(date -u -Iseconds)" "$LOG_FILE" > "$STATUS_FILE"
exec >> "$LOG_FILE" 2>&1

finish() {
    status=$?
    rm -f "$PID_FILE"
    if [ "$status" -eq 0 ]; then
        printf 'SUCCEEDED time=%s\nlog=%s\n' "$(date -u -Iseconds)" "$LOG_FILE" > "$STATUS_FILE"
    else
        printf 'FAILED exit_code=%s time=%s\nlog=%s\n' "$status" "$(date -u -Iseconds)" "$LOG_FILE" > "$STATUS_FILE"
    fi
    exit "$status"
}
trap finish EXIT

set -e
echo "=== AD (ADR+ADF) build started: $(date -u -Iseconds) ==="
echo "log=$LOG_FILE ncomp=$NCOMP source=$SOURCE_DIR"

# wmake shim: bare `-j` (all cores) -> `-j2`
SHIM_DIR=/tmp/wmake-shim
mkdir -p "$SHIM_DIR"
cat > "$SHIM_DIR/wmake" <<'EOF'
#!/usr/bin/env bash
args=()
for a in "$@"; do
    if [ "$a" = "-j" ]; then args+=("-j2"); else args+=("$a"); fi
done
real="${WM_PROJECT_DIR:?WM_PROJECT_DIR unset}/wmake/wmake"
if [ ! -x "$real" ]; then
    IFS=':' read -ra dirs <<< "$PATH"
    for d in "${dirs[@]}"; do
        [ "$d" = /tmp/wmake-shim ] && continue
        if [ -x "$d/wmake" ]; then real="$d/wmake"; break; fi
    done
fi
exec "$real" "${args[@]}"
EOF
chmod +x "$SHIM_DIR/wmake"
export PATH="$SHIM_DIR:$PATH"
hash -r
echo "shim wmake -> $(command -v wmake)"

cd "$SOURCE_DIR"
# NB: OpenFOAM's config.sh/setup breaks under set -e (pop_var_context) -> keep
# errexit OFF while sourcing env scripts, ON while compiling (upstream wrapper
# never enables errexit at all, which is why its failed builds report SUCCEEDED).
# pin_env must re-run after EVERY source: the OpenFOAM bashrc rebuilds PATH and
# drops the wmake shim (verified: `command -v wmake` returns the real wmake
# after sourcing) -> bare `wmake -j` would again use all 10 cores and OOM.
pin_env() {
    export WM_NCOMPPROCS="$NCOMP"
    # move the shim to the VERY FRONT: the OpenFOAM bashrcs prepend their own
    # bin dirs ahead of it, so mere presence in PATH is not enough (PATH lookup
    # then finds the real wmake and -j silently means all 10 cores again)
    local -a dirs=()
    local IFS=':'
    local d
    for d in $PATH; do
        [ "$d" = "$SHIM_DIR" ] || dirs+=("$d")
    done
    unset IFS
    export PATH="$SHIM_DIR:$(IFS=:; printf '%s' "${dirs[*]}")"
    hash -r
}
set +e
source "$DEP_ROOT/loadDAFoam.sh"
set -e
[ -n "$WM_PROJECT" ] || { echo "FATAL: env load failed"; exit 3; }
pin_env

make_mode() {
    echo "=== make_mode WM_AD_MODE=${WM_AD_MODE:-original} WM_NCOMPPROCS=$WM_NCOMPPROCS $(date -u -Iseconds) ==="
    wmakeLnInclude src/adjoint
    ( cd src/adjoint && wmake -j )
    ( cd src/newTurbModels && ./Allmake )
    ( cd src/pyUnitTests && ./Allmake )
    ( cd src/pyDASolvers && ./Allmake )
    if [ -d src/utilities/preProcessing ]; then
        ( cd src/utilities/preProcessing && ./Allmake )
    fi
    echo "=== make_mode done $(date -u -Iseconds) ==="
}

# --- ADR pass ---
sed -i 's/export WM_AD_MODE=.*/export WM_AD_MODE=ADR/' "$DAFOAM_ROOT_PATH/OpenFOAM/OpenFOAM-AD/etc/bashrc"
set +e
source "$DAFOAM_ROOT_PATH/OpenFOAM/OpenFOAM-AD/etc/bashrc"
set -e
[ "$WM_AD_MODE" = "ADR" ] || { echo "FATAL: ADR env load failed"; exit 3; }
pin_env
make_mode

# --- ADF pass ---
set +e
source "$DEP_ROOT/loadDAFoam.sh"
sed -i 's/export WM_AD_MODE=.*/export WM_AD_MODE=ADF/' "$DAFOAM_ROOT_PATH/OpenFOAM/OpenFOAM-AD/etc/bashrc"
source "$DAFOAM_ROOT_PATH/OpenFOAM/OpenFOAM-AD/etc/bashrc"
set -e
[ "$WM_AD_MODE" = "ADF" ] || { echo "FATAL: ADF env load failed"; exit 3; }
pin_env
make_mode

# --- reset to original mode + install python package (ships dafoam/libs/*) ---
set +e
source "$DEP_ROOT/loadDAFoam.sh"
set -e
pin_env
cd "$SOURCE_DIR"
pip install .
ls -R dafoam/libs
echo "=== AD (ADR+ADF) build completed: $(date -u -Iseconds) ==="
echo "*** Build Successful! ***"
