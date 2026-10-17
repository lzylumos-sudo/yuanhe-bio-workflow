#!/usr/bin/env bash
# 用法: ./run.sh <sbvs2 子命令...>   例: ./run.sh doctor
# 解释器优先 $SBVS_PYTHON，其次 PATH 中带 rdkit 的 python3；不写死任何盘符或 /opt/sbvs。
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SBVS_PYTHON:-python3}"
export SBVS_HOME="${SBVS_HOME:-$HERE}" PYTHONPATH="$HERE/source${PYTHONPATH:+:$PYTHONPATH}" PYTHONDONTWRITEBYTECODE=1
exec "$PY" -m sbvs2 "$@"
