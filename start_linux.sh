#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
UI_PYTHON=${YUANHE_UI_PYTHON:-}
if [[ -z "$UI_PYTHON" && -f deployment/runtime.json ]]; then
  UI_PYTHON=$(python3 -c 'import json; print(json.load(open("deployment/runtime.json")).get("chem_python", ""))')
fi
if [[ -z "$UI_PYTHON" ]]; then UI_PYTHON=python3; fi
exec "$UI_PYTHON" -B '00_workbench/run_webui.py' --no-browser --port 8777 "$@"
