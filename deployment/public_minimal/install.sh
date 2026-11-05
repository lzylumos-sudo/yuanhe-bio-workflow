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
lock="$package_dir/linux-64.lock"
if [[ -n ${YUANHE_CONDA_MIRROR:-} ]]; then
  lock="${prefix}.download.lock"
  python3 "$package_dir/../mirror_lock.py" --source "$package_dir/linux-64.lock" \
    --output "$lock" --mirror "$YUANHE_CONDA_MIRROR"
fi
"$mamba" create --no-rc -y -p "$prefix" --file "$lock"
# --isolated ignores PIP_INDEX_URL, so pass the administrator's choice explicitly.
"$prefix/bin/python" -m pip --isolated install --index-url "${PIP_INDEX_URL:-https://pypi.org/simple}" \
  --no-deps --no-build-isolation -r "$package_dir/pip-requirements.txt"
resp_archive=$("$prefix/bin/python" "$package_dir/verify_resp_source.py" --path)
"$prefix/bin/python" -m pip --isolated install --no-deps --no-build-isolation "$resp_archive"
"$prefix/bin/python" -m pip check
echo "Installed independent environment: $prefix"
echo 'Next: run_public.sh ENV_PREFIX NEW_OUTPUT_DIRECTORY'
