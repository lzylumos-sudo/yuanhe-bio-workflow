#!/bin/bash
set -euo pipefail
ROOT="$1"
/opt/sbvs/vendor/micromamba list -p /opt/sbvs/env --explicit >"$ROOT/env/conda-linux-64.lock"
/opt/sbvs/vendor/micromamba env export -p /opt/sbvs/env >"$ROOT/env/environment.yml"
/opt/sbvs/env/bin/python "$ROOT/scripts/freeze_pip.py" "$ROOT"
