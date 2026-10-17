#!/usr/bin/env bash
# Run one named step at a time INSIDE the existing project WSL distribution.
set -Eeuo pipefail
STEP=${1:-help}
if [[ "$STEP" == help ]]; then
  echo 'Usage: bash install_envs.sh base|chem|peptide|reinvent|prep|membrane|gnina|gromacs-cpu|configure|check'
  exit 0
fi
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/../.." && pwd)
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Use Linux x86_64 inside WSL2.'; exit 2; }
[[ -f "$ROOT/02_docking/project.json" ]] || { echo 'Put this bundle in project/deployment/validation_5070.'; exit 2; }
python3 - "$ROOT" "$STEP" <<'PY'
import json,sys,os
from pathlib import Path
r=Path(sys.argv[1])
expected=json.loads((r/'02_docking/project.json').read_text(encoding='utf-8-sig'))['uuid']
markers=[Path('/opt/sbvs/project_uuid'),Path('/opt/yuanhe/project_uuid')]
existing=[m.read_text().strip() for m in markers if m.is_file()]
if any(value!=expected for value in existing):
 raise SystemExit('Conflicting project registration; use the intended WSL distribution.')
if expected not in existing and not (sys.argv[2]=='base' and os.environ.get('YUANHE_ALLOW_NEW_HOST')=='1'):
 raise SystemExit('Unregistered host. For a new dedicated WSL distribution, run base with YUANHE_ALLOW_NEW_HOST=1.')
PY
if [[ "$STEP" == base ]]; then
  sudo apt-get update
  sudo apt-get install -y python3-venv python3-pip git curl ca-certificates build-essential libgomp1 libgl1 libglib2.0-0
  sudo install -d -o "$(id -u)" -g "$(id -g)" /opt/yuanhe
  python3 - "$ROOT" <<'PY'
import json,sys
from pathlib import Path
root=Path(sys.argv[1])
Path('/opt/yuanhe/project_uuid').write_text(json.loads((root/'02_docking/project.json').read_text(encoding='utf-8-sig'))['uuid']+'\n')
PY
  mkdir -p /opt/yuanhe/envs /opt/yuanhe/evidence /opt/yuanhe/cache /opt/yuanhe/python
  [[ -x /opt/yuanhe/tools/bin/python ]] || python3 -m venv /opt/yuanhe/tools
  /opt/yuanhe/tools/bin/python -m pip install 'uv==0.13.1'
fi
[[ -x /opt/yuanhe/tools/bin/uv ]] || { echo 'Run base first.'; exit 2; }
export UV_PYTHON_INSTALL_DIR=/opt/yuanhe/python UV_CACHE_DIR=/opt/yuanhe/cache UV_LINK_MODE=copy
UV=/opt/yuanhe/tools/bin/uv
GEN="$ROOT/01_molecule_generation"
LOG=/opt/yuanhe/evidence/${STEP}_$(date +%Y%m%d_%H%M%S)_$$.log
exec > >(tee -a "$LOG") 2>&1
trap 'echo "FAILED at line $LINENO; keep this log: $LOG"' ERR
newenv() {
  if [[ ! -x "/opt/yuanhe/envs/$1/bin/python" ]]; then
    "$UV" venv --seed --python "$2" "/opt/yuanhe/envs/$1"
  fi
}
freeze() {
  "$UV" pip check --python "/opt/yuanhe/envs/$1/bin/python"
  "$UV" pip freeze --python "/opt/yuanhe/envs/$1/bin/python" > "/opt/yuanhe/evidence/$1.freeze.txt"
}
case "$STEP" in
base) "$UV" --version ;;
chem)
  newenv chem 3.11
  "$UV" pip install --python /opt/yuanhe/envs/chem/bin/python "$ROOT/03_virtual_screening" pytest
  /opt/yuanhe/envs/chem/bin/python -c 'import rdkit,numpy,pandas,scipy,yaml,gemmi,MDAnalysis,prolif,posebusters,zvss; print("Chemical/screening imports OK")'
  freeze chem ;;
peptide)
  newenv peptide2mol 3.10
  P=/opt/yuanhe/envs/peptide2mol/bin/python
  "$P" -m pip install torch==2.7.0 torchvision==0.22.0 torchaudio==2.7.0 --index-url https://download.pytorch.org/whl/cu128
  "$P" -m pip install --no-index --find-links 'https://data.pyg.org/whl/torch-2.7.0+cu128.html' --no-deps \
    'pyg_lib==0.5.0+pt27cu128' 'torch_scatter==2.1.2+pt27cu128' 'torch_sparse==0.6.18+pt27cu128' 'torch_cluster==1.6.3+pt27cu128' 'torch_spline_conv==1.2.2+pt27cu128'
  "$P" -m pip install -r "$HERE/peptide-inference-requirements.txt"
  freeze peptide2mol
  "$P" "$HERE/check_gpu.py" ;;
reinvent)
  # Keep repository lock and CPU inference; its cu126 build is not a 5070 GPU validation.
  UV_PROJECT_ENVIRONMENT=/opt/yuanhe/envs/reinvent4 "$UV" sync --project "$GEN/vendor/REINVENT4" --python 3.11 --frozen --no-editable
  /opt/yuanhe/envs/reinvent4/bin/reinvent --help
  freeze reinvent4 ;;
prep)
  newenv prep 3.11
  "$UV" pip install --python /opt/yuanhe/envs/prep/bin/python 'pdbfixer==1.12.0' openmm
  /opt/yuanhe/envs/prep/bin/python -c 'import pdbfixer,openmm; print(openmm.__version__)'
  freeze prep ;;
membrane) bash "$ROOT/deployment/install_membrane.sh" ;;
gnina)
  source "$ROOT/02_docking/scripts/runtime_env.sh"
  [[ -x /opt/sbvs/vendor/gnina ]] || { echo 'Existing project runtime is missing GNINA; registration/repair required.'; exit 3; }
  /opt/sbvs/vendor/gnina --version
  sha256sum /opt/sbvs/vendor/gnina
  echo 'Version check only. GPU CNN execution must still be tested on target.' ;;
gromacs-cpu)
  sudo apt-get update
  sudo apt-get install -y gromacs
  /usr/bin/gmx --version
  echo 'Distribution GROMACS installed. Inspect GPU support; do not label a CPU build CUDA-enabled.' ;;
configure)
  python3 - "$ROOT" <<'PY'
import json,sys,shutil,time
from pathlib import Path
r=Path(sys.argv[1]);p=r/'deployment/runtime.json'
cfg=json.loads(p.read_text()) if p.exists() else {}
if p.exists():shutil.copyfile(p,p.with_name('runtime.before-validation-'+str(time.time_ns())+'.json'))
cfg.update(chem_python='/opt/yuanhe/envs/chem/bin/python',screen_python='/opt/yuanhe/envs/chem/bin/python',peptide_python='/opt/yuanhe/envs/peptide2mol/bin/python',reinvent_executable='/opt/yuanhe/envs/reinvent4/bin/reinvent',prep_python='/opt/yuanhe/envs/prep/bin/python',md_python='/opt/yuanhe/envs/membrane/bin/python',gnina='/opt/sbvs/vendor/gnina',gmx=cfg.get('gmx') or '/usr/bin/gmx')
p.write_text(json.dumps(cfg,ensure_ascii=False,indent=2))
print('Configured runtime paths. Existing runs are not changed.')
PY
  ;;
check)
  python3 "$ROOT/deployment/doctor.py" --config "$ROOT/deployment/runtime.json" --out /opt/yuanhe/evidence/doctor.json
  /opt/yuanhe/envs/peptide2mol/bin/python "$HERE/check_gpu.py" ;;
*) echo "Unknown step: $STEP"; exit 2 ;;
esac
echo "STEP COMPLETE: $STEP. Log: $LOG. No target inference/docking/MD was launched."
