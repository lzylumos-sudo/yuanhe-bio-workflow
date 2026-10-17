#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/runtime.sh"
[[ -x "$NEXUS_ENV_PREFIX/bin/python" ]] || { echo 'Environment missing: run bash setup.sh'; exit 2; }
exec "$NEXUS_ENV_PREFIX/bin/python" "$PROJECT_ROOT/scripts/run_md.py" "$@"
