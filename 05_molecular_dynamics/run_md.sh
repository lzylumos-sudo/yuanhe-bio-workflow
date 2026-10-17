#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  echo 'Execution target: x86-64 Ubuntu/Linux. Mac may run on Python for input-layer tests; this does not establish Linux MD validation.' >&2
  exit 2
fi
source "$PROJECT_ROOT/runtime.sh"
[[ -x "$NEXUS_ENV_PREFIX/bin/python" ]] || { echo 'First run in Ubuntu using bash setup.sh'; exit 2; }
exec "$NEXUS_ENV_PREFIX/bin/python" "$PROJECT_ROOT/scripts/run_workflow.py" "$@"
