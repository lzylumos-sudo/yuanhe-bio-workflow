#!/usr/bin/env bash
# Optional portable CUDA build; never installs drivers or changes a site service.
set -Eeuo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo 'Usage: bash deployment/build_gromacs.sh /absolute/NEW/prefix /absolute/cuda [jobs]' >&2
  exit 2
fi
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Linux x86_64 required.' >&2; exit 2; }
BUILD_PREFIX=$1
CUDA_ROOT=$2
BUILD_JOBS=${3:-4}
[[ "$BUILD_PREFIX" == /* && ! -e "$BUILD_PREFIX" ]] || { echo 'Choose a new absolute prefix; existing files are preserved.' >&2; exit 2; }
[[ "$CUDA_ROOT" == /* && -x "$CUDA_ROOT/bin/nvcc" ]] || { echo 'Provide an installed CUDA toolkit, not only a driver.' >&2; exit 2; }
[[ "$BUILD_JOBS" =~ ^[0-9]+$ && "$BUILD_JOBS" -ge 1 && "$BUILD_JOBS" -le 128 ]] || exit 2
BUILD_CMAKE=${YUANHE_CMAKE:-cmake}
BUILD_CC=$(command -v "${YUANHE_CC:-gcc}")
BUILD_CXX=$(command -v "${YUANHE_CXX:-g++}")
for tool in "$BUILD_CMAKE" curl tar sha256sum make; do command -v "$tool" >/dev/null; done
mkdir -p "$BUILD_PREFIX/evidence"
exec > >(tee -a "$BUILD_PREFIX/evidence/build.log") 2>&1
trap 'echo "Build failed at line $LINENO; keep logs and use a new prefix after resolving the cause."' ERR
"$BUILD_CMAKE" --version
"$CUDA_ROOT/bin/nvcc" --version
"$BUILD_CXX" --version
SOURCE_ARCHIVE="$BUILD_PREFIX/gromacs-2026.4.tar.gz"
curl -fL --retry 3 --connect-timeout 20 --max-time 1800 \
  -o "$SOURCE_ARCHIVE" https://ftp.gromacs.org/gromacs/gromacs-2026.4.tar.gz
printf '%s  %s\n' '58eda60979b124fcf9dff8b1b63cedfbc72fdb36d99e921dc1bdc310571098b6' "$SOURCE_ARCHIVE" | sha256sum -c -
tar --no-same-owner -xzf "$SOURCE_ARCHIVE" -C "$BUILD_PREFIX"
"$BUILD_CMAKE" -S "$BUILD_PREFIX/gromacs-2026.4" -B "$BUILD_PREFIX/build" \
  -DCMAKE_INSTALL_PREFIX="$BUILD_PREFIX/install" \
  -DCMAKE_C_COMPILER="$BUILD_CC" -DCMAKE_CXX_COMPILER="$BUILD_CXX" \
  -DCMAKE_CUDA_COMPILER="$CUDA_ROOT/bin/nvcc" -DCMAKE_CUDA_HOST_COMPILER="$BUILD_CXX" \
  -DCUDAToolkit_ROOT="$CUDA_ROOT" -DGMX_GPU=CUDA -DGMX_MPI=OFF \
  -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
"$BUILD_CMAKE" --build "$BUILD_PREFIX/build" --parallel "$BUILD_JOBS"
"$BUILD_CMAKE" --build "$BUILD_PREFIX/build" --target check --parallel "$BUILD_JOBS"
"$BUILD_CMAKE" --install "$BUILD_PREFIX/build"
"$BUILD_PREFIX/install/bin/gmx" --version | tee "$BUILD_PREFIX/evidence/gmx-version.txt"
python3 - "$BUILD_PREFIX/evidence/gmx-version.txt" <<'PY'
from pathlib import Path
import sys
assert any('GPU support:' in line and 'CUDA' in line for line in Path(sys.argv[1]).read_text().splitlines())
PY
echo "CUDA build and upstream checks completed. Use --gmx $BUILD_PREFIX/install/bin/gmx in deploy_workbench.py."
echo 'Target membrane dynamics and checkpoint recovery still require real acceptance.'
