#!/usr/bin/env bash
# Administrator-only installation into a dedicated prefix. No driver/system changes.
set -Eeuo pipefail
if [[ $# -ne 2 ]]; then
  echo 'Usage: bash deployment/install_workbench.sh runtime|peptide|reinvent|vina|configure|check /absolute/prefix' >&2
  exit 2
fi
STEP=$1
PREFIX=$2
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Requires Linux x86_64 / WSL2.' >&2; exit 2; }
[[ "$PREFIX" == /* && "$PREFIX" != / ]] || { echo 'Use an absolute dedicated prefix.' >&2; exit 2; }
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/.." && pwd)
mkdir -p "$PREFIX/evidence"
LOG="$PREFIX/evidence/${STEP}_$(date +%Y%m%d_%H%M%S)_$$.log"
exec > >(tee -a "$LOG") 2>&1
trap 'echo "FAILED at line $LINENO. Preserve $LOG; do not report deployment complete."' ERR
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
RUNTIME="$PREFIX/envs/runtime"
TOOLS="$PREFIX/tools"
UV="$TOOLS/bin/uv"
GEN="$ROOT/01_molecule_generation"
export UV_PYTHON_INSTALL_DIR="$PREFIX/python" UV_CACHE_DIR="$PREFIX/cache" UV_LINK_MODE=copy
freeze() {
  if [[ -x "$UV" ]]; then
    "$UV" pip check --python "$1/bin/python"
    "$UV" pip freeze --python "$1/bin/python" > "$PREFIX/evidence/$2.freeze.txt"
  else
    "$1/bin/python" -m pip check
    "$1/bin/python" -m pip freeze > "$PREFIX/evidence/$2.freeze.txt"
  fi
}
case "$STEP" in
runtime)
  # Reuse the independently reconstructed Linux lock; do not merge unrelated envs.
  bash "$HERE/public_minimal/install.sh" "$RUNTIME"
  "$RUNTIME/bin/python" -m pip install --no-deps --no-build-isolation "$ROOT/03_virtual_screening"
  freeze "$RUNTIME" runtime
  "$RUNTIME/bin/python" -m venv "$TOOLS"
  "$TOOLS/bin/python" -m pip install 'uv==0.13.1'
  cp "$HERE/public_minimal/linux-64.lock" "$PREFIX/evidence/runtime-linux-64.lock"
  cp "$HERE/public_minimal/pip-requirements.txt" "$PREFIX/evidence/runtime-pip-requirements.txt"
  ;;
peptide)
  [[ -x "$UV" ]] || { echo 'Install runtime first.' >&2; exit 2; }
  [[ -d "$GEN/vendor/Peptide2Mol" && -d "$GEN/vendor/REINVENT4" ]] || { echo 'Restore both pinned vendor sources and weights first; see WORKBENCH_DEPLOYMENT.md.' >&2; exit 2; }
  "$RUNTIME/bin/python" "$HERE/apply_vendor_patches.py"
  P="$PREFIX/envs/peptide2mol"
  [[ ! -e "$P" ]] || { echo 'Environment exists. Inspect the retained log; use a new prefix after a failed install.' >&2; exit 2; }
  "$UV" venv --seed --python 3.10 "$P"
  "$P/bin/python" -m pip install torch==2.7.0 torchvision==0.22.0 torchaudio==2.7.0 --index-url https://download.pytorch.org/whl/cu128
  "$P/bin/python" -m pip install --no-index --find-links 'https://data.pyg.org/whl/torch-2.7.0+cu128.html' --no-deps \
    'pyg_lib==0.5.0+pt27cu128' 'torch_scatter==2.1.2+pt27cu128' 'torch_sparse==0.6.18+pt27cu128' 'torch_cluster==1.6.3+pt27cu128' 'torch_spline_conv==1.2.2+pt27cu128'
  "$P/bin/python" -m pip install -r "$HERE/validation_5070/peptide-inference-requirements.txt"
  freeze "$P" peptide2mol
  "$P/bin/python" "$HERE/validation_5070/check_gpu.py"
  ;;
reinvent)
  [[ -x "$UV" && -f "$GEN/vendor/REINVENT4/uv.lock" ]] || { echo 'Install runtime and restore pinned REINVENT4 including uv.lock first.' >&2; exit 2; }
  "$RUNTIME/bin/python" "$HERE/apply_vendor_patches.py"
  [[ ! -e "$PREFIX/envs/reinvent4" ]] || { echo 'Environment exists; refusing an implicit upgrade.' >&2; exit 2; }
  UV_PROJECT_ENVIRONMENT="$PREFIX/envs/reinvent4" "$UV" sync --project "$GEN/vendor/REINVENT4" --python 3.11 --frozen
  "$PREFIX/envs/reinvent4/bin/reinvent" --help
  freeze "$PREFIX/envs/reinvent4" reinvent4
  ;;
vina)
  SBVS_PROJECT_ROOT="$ROOT/02_docking" YUANHE_VENDOR_DIR="$PREFIX/bin" "$RUNTIME/bin/python" "$ROOT/02_docking/scripts/install_vendor.py" vina
  ;;
configure)
  [[ -x "$RUNTIME/bin/python" ]] || { echo 'Install runtime first.' >&2; exit 2; }
  [[ -n ${YUANHE_GMX:-} && -x "$YUANHE_GMX" ]] || { echo 'Set YUANHE_GMX to the absolute CUDA-enabled GROMACS executable.' >&2; exit 2; }
  "$YUANHE_GMX" --version
  "$RUNTIME/bin/python" - "$ROOT" "$PREFIX" "$YUANHE_GMX" <<'PY'
import json,shutil,sys,time
from pathlib import Path
root,prefix=map(Path,sys.argv[1:3]);dest=root/'deployment/runtime.json'
cfg=json.loads(dest.read_text()) if dest.is_file() else {}
if dest.exists():shutil.copy2(dest,prefix/'evidence'/('runtime.before-'+str(time.time_ns())+'.json'))
p=str(prefix/'envs/runtime/bin/python')
cfg.update(chem_python=p,screen_python=p,prep_python=p,md_python=p,
           peptide_python=str(prefix/'envs/peptide2mol/bin/python'),
           reinvent_executable=str(prefix/'envs/reinvent4/bin/reinvent'),
           vina=str(prefix/'bin/vina'),gmx=sys.argv[3])
# Explicit research defaults let administrators configure evidence before users create runs.
if 'research' not in cfg:
    cfg['research']=json.loads((root/'00_workbench/configs/research_dual_blackwell.json').read_text())
dest.write_text(json.dumps(cfg,indent=2)+'\n')
print('Runtime paths saved. Configure research.predocking evidence, then run check. No task started.')
PY
  ;;
check)
  "$RUNTIME/bin/python" "$HERE/doctor.py" --config "$HERE/runtime.json" --out "$PREFIX/evidence/readiness.json"
  ;;
*) echo "Unknown step: $STEP" >&2; exit 2 ;;
esac
echo "Step completed: $STEP. Installation logs: $LOG. Full research acceptance remains separate."
