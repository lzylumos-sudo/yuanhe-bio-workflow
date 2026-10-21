#!/usr/bin/env bash
# Creates a NEW independent environment. Never activates or modifies base.
set -euo pipefail
package_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ $# -ne 1 ]]; then echo "Usage: MICROMAMBA=/path/to/micromamba bash install.sh /absolute/new/env" >&2; exit 2; fi
prefix="$1"
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || { echo 'Linux x86_64 / WSL2 required' >&2; exit 2; }
[[ "$prefix" == /* && ! -e "$prefix" ]] || { echo 'Environment path must be absolute and must not exist' >&2; exit 2; }
mamba="${MICROMAMBA:-micromamba}"
command -v "$mamba" >/dev/null || { echo 'Install micromamba first and set MICROMAMBA to its executable' >&2; exit 2; }
export MAMBA_ROOT_PREFIX="${prefix}-mamba-root"
export CONDA_PKGS_DIRS="${YUANHE_PACKAGE_CACHE:-${prefix}-packages}"
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
"$mamba" create --no-rc -y -p "$prefix" --file "$package_dir/linux-64.lock"
"$prefix/bin/python" -m pip --isolated install --no-deps --no-build-isolation -r "$package_dir/pip-requirements.txt"
"$prefix/bin/python" -m pip check
echo "Installed independent environment: $prefix"
echo 'Next: run_public.sh ENV_PREFIX NEW_OUTPUT_DIRECTORY'
