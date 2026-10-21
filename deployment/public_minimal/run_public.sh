#!/usr/bin/env bash
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ $# -ne 2 ]]; then echo 'Usage: bash run_public.sh ENV_PREFIX NEW_OUTPUT_DIRECTORY' >&2; exit 2; fi
prefix="$1"
output="$2"
[[ "$prefix" == /* && -x "$prefix/bin/python" ]] || { echo 'Absolute installed environment required' >&2; exit 2; }
[[ "$output" == /* && ! -e "$output" ]] || { echo 'Output path must be absolute and NEW' >&2; exit 2; }
export PATH="$prefix/bin:/usr/local/bin:/usr/bin:/bin"
export CONDA_PREFIX="$prefix" AMBERHOME="$prefix" PYTHONDONTWRITEBYTECODE=1
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENMM_CPU_THREADS=1 MPLBACKEND=Agg
mkdir -p -- "$output"
export MPLCONFIGDIR="$output/cache/matplotlib" XDG_CACHE_HOME="$output/cache" PSI_SCRATCH="$output/scratch"
mkdir -p -- "$PSI_SCRATCH"
"$prefix/bin/python" -u "$package_dir/run_public.py" "$output"
