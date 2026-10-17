#!/bin/bash
set -euo pipefail
export SBVS_PROJECT_ROOT="${SBVS_PROJECT_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
if [ ! -x /opt/sbvs/env/bin/python ] || [ ! -f /opt/sbvs/project_uuid ]; then
 echo 'Use START_WINDOWS.ps1 or SBVS.ps1 in the project dedicated WSL. This distribution is not configured.' >&2; exit 78
fi
export PYTHONPATH="$SBVS_PROJECT_ROOT/src"
export PLAYWRIGHT_BROWSERS_PATH=/opt/sbvs/browsers
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source "$SBVS_PROJECT_ROOT/scripts/runtime_env.sh"
exec /opt/sbvs/env/bin/python -m sbvs.cli "$@"
