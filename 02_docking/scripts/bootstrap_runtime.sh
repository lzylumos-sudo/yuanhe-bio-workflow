#!/bin/bash
set -euo pipefail
ROOT="$1"
UUID="$2"
export DEBIAN_FRONTEND=noninteractive
mkdir -p /opt/sbvs/{vendor,cache,work,state,browsers}
printf '%s' "$UUID" >/opt/sbvs/project_uuid
apt-get update
apt-get install -y --no-install-recommends ca-certificates curl bzip2 unzip libgomp1 libglib2.0-0 libgl1 libxrender1 libxext6 procps
id sbvs >/dev/null 2>&1 || useradd -m -s /bin/bash sbvs
printf '[user]\ndefault=sbvs\n[automount]\nenabled=true\noptions="metadata"\n' >/etc/wsl.conf
chown -R sbvs:sbvs /opt/sbvs
if [ ! -x /opt/sbvs/vendor/micromamba ]; then
 curl -fL --retry 2 --connect-timeout 20 --max-time 300 https://micro.mamba.pm/api/micromamba/linux-64/2.3.3 -o /opt/sbvs/cache/micromamba.tar.bz2
 tar -xjf /opt/sbvs/cache/micromamba.tar.bz2 -C /opt/sbvs/vendor --strip-components=1 bin/micromamba
fi
export MAMBA_ROOT_PREFIX=/opt/sbvs/mamba
if [ ! -x /opt/sbvs/env/bin/python ]; then
 if [ -s "$ROOT/env/conda-linux-64.lock" ]; then
  /opt/sbvs/vendor/micromamba create -y -p /opt/sbvs/env -f "$ROOT/env/conda-linux-64.lock"
 else
  /opt/sbvs/vendor/micromamba create -y -p /opt/sbvs/env -c conda-forge --strict-channel-priority python=3.11 rdkit meeko=0.8.0 gemmi biopython numpy scipy pandas pyyaml flask waitress pytest psutil filelock requests pip
 fi
fi
export PIP_CACHE_DIR=/opt/sbvs/cache/pip PLAYWRIGHT_BROWSERS_PATH=/opt/sbvs/browsers
/opt/sbvs/env/bin/python -m pip install playwright==1.55.0 molscrub==0.2.2
/opt/sbvs/env/bin/python -m playwright install --with-deps chromium
bash "$ROOT/scripts/refresh_locks.sh" "$ROOT"
chown -R sbvs:sbvs /opt/sbvs
/opt/sbvs/env/bin/python -c 'import rdkit,meeko,Bio; print(rdkit.__version__,meeko.__version__,Bio.__version__)'
