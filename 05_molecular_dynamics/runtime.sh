#!/usr/bin/env bash
# Source this file; no environment is shared with another module.
PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
NEXUS_CACHE="${XDG_CACHE_HOME:-$HOME/.cache}/nexus-structural-biology-v1"
export NEXUS_ENV_PREFIX="${NEXUS_ENV_PREFIX:-$NEXUS_CACHE/environments/md}"
export AMBERHOME="$NEXUS_ENV_PREFIX"
export PATH="$NEXUS_ENV_PREFIX/bin:$PATH"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENMM_CPU_THREADS="${OPENMM_CPU_THREADS:-4}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$NEXUS_CACHE/cache/matplotlib}"
export PYTHONPYCACHEPREFIX="${PYTHONPYCACHEPREFIX:-$NEXUS_CACHE/cache/pycache}"
mkdir -p "$MPLCONFIGDIR" "$PYTHONPYCACHEPREFIX"
