#!/usr/bin/env bash
# ZLIS launcher (macOS / Linux). Works from any folder, including paths with spaces or non-ASCII characters.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${ZLIS_ENV:-zlis}"

find_python() {
  for root in "${CONDA_PREFIX_1:-}" "${CONDA_PREFIX:-}" "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/mambaforge" \
              "/opt/miniconda3" "/opt/homebrew/Caskroom/miniforge/base" "/usr/local/miniconda3"; do
    [ -n "$root" ] && [ -x "$root/envs/$ENV_NAME/bin/python3" ] && { echo "$root/envs/$ENV_NAME/bin/python3"; return 0; }
  done
  if command -v conda >/dev/null 2>&1; then
    base="$(conda info --base 2>/dev/null)"
    [ -x "$base/envs/$ENV_NAME/bin/python3" ] && { echo "$base/envs/$ENV_NAME/bin/python3"; return 0; }
  fi
  return 1
}

PY="$(find_python)" || {
  echo "ZLIS environment '$ENV_NAME' not found. Run: bash \"$HERE/env/setup_mac.sh\"" >&2
  exit 2
}
export PYTHONPATH="$HERE/src${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONWARNINGS=ignore
export PYTHONDONTWRITEBYTECODE=1
exec "$PY" -m zlis "$@"
