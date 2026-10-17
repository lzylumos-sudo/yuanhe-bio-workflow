#!/usr/bin/env bash
# Ubuntu / WSL2 Ubuntu installer. Existing runs and runtime archives are untouched.
set -euo pipefail
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  echo 'This runtime targets x86-64 Ubuntu/Linux (including WSL2 Ubuntu).'
  echo 'Reading or inspecting the project on a Mac is supported; native MD installation must run on Linux.'
  exit 2
fi
repair=false
if [[ "${1:-}" == --repair && $# -eq 1 ]]; then
  repair=true
elif [[ $# -ne 0 ]]; then
  echo 'Usage: bash setup.sh [--repair]' >&2
  exit 2
fi
source "$(dirname -- "${BASH_SOURCE[0]}")/runtime.sh"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
evidence="$PROJECT_ROOT/installation_evidence"
mkdir -p "$evidence" "$(dirname -- "$NEXUS_ENV_PREFIX")"
lock_dir="${NEXUS_ENV_PREFIX}.install.lock"
if ! mkdir "$lock_dir" 2>/dev/null; then
  echo "Another installation may be active: $lock_dir" >&2
  echo 'Check for a running installer before removing a stale empty installation lock.' >&2
  exit 2
fi
trap 'rmdir -- "$lock_dir" 2>/dev/null || true' EXIT
check_environment() {
  local report="$evidence/environment-check-${stamp}.json"
  if [[ ! -x "$NEXUS_ENV_PREFIX/bin/python" ]]; then
    return 1
  fi
  "$NEXUS_ENV_PREFIX/bin/python" "$PROJECT_ROOT/scripts/check_environment.py" --native --output "$report"
}
if [[ -e "$NEXUS_ENV_PREFIX" ]]; then
  if check_environment; then
    echo 'Existing MD environment passed dependency, executable and native CPU/checkpoint checks.'
    echo 'Next run the v2 molecular smoke test; installation checks do not validate ligand parameterization.'
    exit 0
  fi
  if [[ "$repair" != true ]]; then
    echo "Existing environment is incomplete or failed validation: $NEXUS_ENV_PREFIX" >&2
    echo 'Run bash setup.sh --repair to preserve it under a timestamped name and restore a new environment.' >&2
    exit 2
  fi
  backup="${NEXUS_ENV_PREFIX}.before-repair-${stamp}-$$"
  mv -- "$NEXUS_ENV_PREFIX" "$backup"
  echo "Previous environment preserved at: $backup"
fi
archive="$PROJECT_ROOT/runtime/environment-linux-64.tar.gz"
if [[ -f "$archive" ]]; then
  checksum="$archive.sha256"
  [[ -f "$checksum" ]] || { echo 'Archive checksum missing; restore refused.' >&2; exit 2; }
  (cd "$PROJECT_ROOT/runtime" && sha256sum -c environment-linux-64.tar.gz.sha256)
  stage="${NEXUS_ENV_PREFIX}.extract-${stamp}-$$"
  mkdir "$stage"
  echo "Extracting the Linux environment into: $stage"
  tar -xzf "$archive" -C "$stage"
  [[ -x "$stage/bin/python" && -f "$stage/bin/conda-unpack" ]] || {
    echo "Extracted archive is incomplete; preserved for inspection at $stage" >&2; exit 2;
  }
  mv -- "$stage" "$NEXUS_ENV_PREFIX"
  # conda-unpack must run after the environment reaches its final prefix.
  "$NEXUS_ENV_PREFIX/bin/python" "$NEXUS_ENV_PREFIX/bin/conda-unpack"
else
  echo 'No bundled archive found; rebuilding from pinned explicit package URLs requires Internet access.'
  lockfile="$PROJECT_ROOT/env/linux-64-explicit.lock"
  [[ -f "$lockfile" ]] || { echo 'Pinned explicit lock missing; restore the lock or bundled archive.' >&2; exit 2; }
  # Refuse the malformed legacy environment-list header or unrelated content.
  awk '
    /^[[:space:]]*$/ || /^#/ {next}
    !seen {if ($0 != "@EXPLICIT") exit 1; seen=1; next}
    !/^https:\/\/conda\.anaconda\.org\/conda-forge\/(linux-64|noarch)\/[^[:space:]]+\.(conda|tar\.bz2)$/ {exit 1}
    END {if (!seen) exit 1}
  ' "$lockfile" || { echo 'Explicit package lock is malformed; rebuild refused.' >&2; exit 2; }
  if ! command -v micromamba >/dev/null; then
    command -v curl >/dev/null || { echo 'Install curl first (sudo apt install curl).' >&2; exit 2; }
    mkdir -p "$NEXUS_CACHE/bootstrap"
    curl --fail --location --retry 3 https://micro.mamba.pm/api/micromamba/linux-64/latest -o "$NEXUS_CACHE/bootstrap/micromamba.tar.bz2"
    tar -xjf "$NEXUS_CACHE/bootstrap/micromamba.tar.bz2" -C "$NEXUS_CACHE/bootstrap" bin/micromamba
    export PATH="$NEXUS_CACHE/bootstrap/bin:$PATH"
  fi
  # A conda-created environment is installed directly at its permanent prefix;
  # moving it afterward would invalidate paths embedded in package files.
  micromamba create --yes --prefix "$NEXUS_ENV_PREFIX" --file "$lockfile"
fi
if ! check_environment; then
  echo 'Installation validation failed. Diagnostic JSON is in installation_evidence.' >&2
  echo 'The partial environment is preserved; inspect the error, then use bash setup.sh --repair if restoration is needed.' >&2
  exit 2
fi
echo 'MD environment passed full imports, executable checks and native CPU/checkpoint smoke.'
echo 'Next run the v2 molecular smoke test; no production simulation has been started.'
