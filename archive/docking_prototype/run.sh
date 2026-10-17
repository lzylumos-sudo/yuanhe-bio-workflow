#!/usr/bin/env bash
# Usage: ./run.sh <sbvs2 Subcommands...>   Example: ./run.sh doctor
# Interpreter preference: $SBVS_PYTHON, then PATH containing rdkit associated python3; No hard-coded drive letter or /opt/sbvs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SBVS_PYTHON:-python3}"
export SBVS_HOME="${SBVS_HOME:-$HERE}" PYTHONPATH="$HERE/source${PYTHONPATH:+:$PYTHONPATH}" PYTHONDONTWRITEBYTECODE=1
exec "$PY" -m sbvs2 "$@"
