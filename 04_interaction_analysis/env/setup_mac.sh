#!/usr/bin/env bash
# Creates the 'zlis' conda environment on this Mac. The environment lives on the
# internal disk (conda default), never on an exFAT drive.
set -e
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if ! command -v conda >/dev/null 2>&1; then
  for c in "$HOME/miniconda3/bin/conda" "/opt/miniconda3/bin/conda" "$HOME/miniforge3/bin/conda"; do
    [ -x "$c" ] && { export PATH="$(dirname "$c"):$PATH"; break; }
  done
fi
command -v conda >/dev/null 2>&1 || { echo "conda not found. Install Miniforge first: https://github.com/conda-forge/miniforge" >&2; exit 2; }
conda env create -f "$HERE/environment.yml" --yes || conda env update -n zlis -f "$HERE/environment.yml" --prune
echo "Done. Start with: bash run.sh web"
