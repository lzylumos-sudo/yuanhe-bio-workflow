#!/usr/bin/env bash
# Optional: compile a separate CUDA-enabled GROMACS in the project WSL distro.
# Requires the NVIDIA WSL CUDA 12.8 TOOLKIT already installed (never Linux drivers).
set -Eeuo pipefail
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || exit 2
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/../.." && pwd)
python3 - "$ROOT" <<'PYREG'
import json,sys
from pathlib import Path
expected=json.loads((Path(sys.argv[1])/'02_docking/project.json').read_text(encoding='utf-8-sig'))['uuid']
markers=[Path('/opt/sbvs/project_uuid'),Path('/opt/yuanhe/project_uuid')]
values=[p.read_text().strip() for p in markers if p.is_file()]
if expected not in values or any(v!=expected for v in values):
    raise SystemExit('Use the registered dedicated project host; run the base installation first.')
PYREG
[[ -x /usr/local/cuda-12.8/bin/nvcc ]] || { echo 'Install CUDA Toolkit 12.8 for WSL first; see guide.'; exit 2; }
[[ -x /opt/yuanhe/tools/bin/python ]] || { echo 'Run install_envs.sh base first.'; exit 2; }
sudo apt-get update
sudo apt-get install -y gcc-12 g++-12 build-essential libfftw3-dev
/opt/yuanhe/tools/bin/python -m pip install 'cmake==3.31.6'
mkdir -p /opt/yuanhe/build /opt/yuanhe/evidence
cd /opt/yuanhe/build
exec > >(tee -a "/opt/yuanhe/evidence/gromacs_build_$(date +%Y%m%d_%H%M%S)_$$.log") 2>&1
curl -fL --retry 3 -o gromacs-2026.4.tar.gz https://ftp.gromacs.org/gromacs/gromacs-2026.4.tar.gz
echo 'fbe7d5071dee17e3f448502b1fd30d55  gromacs-2026.4.tar.gz' | md5sum -c -
sha256sum gromacs-2026.4.tar.gz > /opt/yuanhe/evidence/gromacs-source.sha256
[[ -d gromacs-2026.4 ]] || tar -xzf gromacs-2026.4.tar.gz
CMAKE=/opt/yuanhe/tools/bin/cmake
"$CMAKE" -S gromacs-2026.4 -B gromacs-2026.4-build \
 -DCMAKE_INSTALL_PREFIX=/opt/yuanhe/gromacs-2026.4 \
 -DCMAKE_C_COMPILER=/usr/bin/gcc-12 -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 \
 -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.8/bin/nvcc \
 -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-12 -DCUDAToolkit_ROOT=/usr/local/cuda-12.8 \
 -DGMX_GPU=CUDA -DGMX_SIMD=SSE2 -DGMX_MPI=OFF -DGMX_BUILD_OWN_FFTW=OFF \
 -DREGRESSIONTEST_DOWNLOAD=ON
"$CMAKE" --build gromacs-2026.4-build --parallel 4
"$CMAKE" --build gromacs-2026.4-build --target check --parallel 2
"$CMAKE" --install gromacs-2026.4-build
/opt/yuanhe/gromacs-2026.4/bin/gmx --version
echo 'Build completed. Configure gmx=/opt/yuanhe/gromacs-2026.4/bin/gmx; target membrane MD is not yet validated.'
