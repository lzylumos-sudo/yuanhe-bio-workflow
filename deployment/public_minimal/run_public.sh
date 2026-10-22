#!/usr/bin/env bash
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
resume=()
if [[ "${1:-}" == "--resume" ]]; then resume=(--resume); shift; fi
if [[ $# -ne 2 ]]; then echo 'Usage: bash run_public.sh [--resume] ENV_PREFIX OUTPUT_DIRECTORY' >&2; exit 2; fi
prefix="$1"
output="$2"
[[ "$prefix" == /* && -x "$prefix/bin/python" ]] || { echo 'Absolute installed environment required' >&2; exit 2; }
[[ "$output" == /* ]] || { echo 'Output path must be absolute' >&2; exit 2; }
if [[ ${#resume[@]} -eq 0 ]]; then
  [[ ! -e "$output" ]] || { echo 'Output exists; use --resume for a recorded run' >&2; exit 2; }
  mkdir -- "$output"
else
  [[ -d "$output" && -f "$output/RUN_STATE.json" ]] || { echo 'Recorded output with RUN_STATE.json required' >&2; exit 2; }
fi
export PATH="$prefix/bin:/usr/local/bin:/usr/bin:/bin"
export CONDA_PREFIX="$prefix" AMBERHOME="$prefix" PYTHONDONTWRITEBYTECODE=1
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENMM_CPU_THREADS=1 MPLBACKEND=Agg
export MPLCONFIGDIR="$output/cache/matplotlib" XDG_CACHE_HOME="$output/cache" PSI_SCRATCH="$output/scratch"
mkdir -p -- "$PSI_SCRATCH"
exec "$prefix/bin/python" -u "$package_dir/run_public.py" "$output" "${resume[@]}"
