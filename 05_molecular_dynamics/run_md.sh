#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  echo '运行目标为 x86-64 Ubuntu/Linux。Mac 上可用 Python 执行输入层测试，不代表 Linux MD 已验证。' >&2
  exit 2
fi
source "$PROJECT_ROOT/runtime.sh"
[[ -x "$NEXUS_ENV_PREFIX/bin/python" ]] || { echo '请先在 Ubuntu 中运行 bash setup.sh'; exit 2; }
exec "$NEXUS_ENV_PREFIX/bin/python" "$PROJECT_ROOT/scripts/run_workflow.py" "$@"
