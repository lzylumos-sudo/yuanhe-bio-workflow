#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if [[ -d /opt/sbvs/env && -f 02_docking/scripts/runtime_env.sh ]]; then
  source 02_docking/scripts/runtime_env.sh
fi
exec "${YUANHE_UI_PYTHON:-python3}" -B '00_workbench/run_webui.py' --no-browser --port 8777
