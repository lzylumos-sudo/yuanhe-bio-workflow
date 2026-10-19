#!/bin/bash
set -euo pipefail
ROOT="$1"
SBVS_PYTHON="${SBVS_ENV_PREFIX:-/opt/sbvs/env}/bin/python"
"$SBVS_PYTHON" "$ROOT/scripts/export_environment.py" "$ROOT"
"$SBVS_PYTHON" "$ROOT/scripts/freeze_pip.py" "$ROOT"
