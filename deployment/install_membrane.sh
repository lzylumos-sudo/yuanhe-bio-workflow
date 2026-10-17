#!/usr/bin/env bash
# Dedicated environment only. Never modifies the chemistry or model environments.
set -Eeuo pipefail
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Requires Linux x86_64 / WSL2.'; exit 2; }
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
[[ -d /opt/yuanhe && -w /opt/yuanhe ]] || { echo 'Complete the base installation first; /opt/yuanhe must be writable by your WSL user.'; exit 2; }
mkdir -p /opt/yuanhe/bin /opt/yuanhe/cache /opt/yuanhe/evidence
export MAMBA_ROOT_PREFIX=/opt/yuanhe/mamba
MM=/opt/yuanhe/bin/micromamba
if [[ ! -x "$MM" ]]; then
  archive=/opt/yuanhe/cache/micromamba.tar.bz2
  curl --fail --location --retry 3 'https://micro.mamba.pm/api/micromamba/linux-64/latest' --output "$archive"
  sha256sum "$archive" > /opt/yuanhe/evidence/micromamba-download.sha256
  tar -xjf "$archive" -C /opt/yuanhe bin/micromamba
fi
"$MM" --version > /opt/yuanhe/evidence/micromamba.version.txt
ENV=/opt/yuanhe/envs/membrane
if [[ ! -x "$ENV/bin/python" ]]; then
  if [[ -f /opt/yuanhe/evidence/membrane-linux-64.explicit.txt ]]; then
    "$MM" create --yes --no-rc --prefix "$ENV" --file /opt/yuanhe/evidence/membrane-linux-64.explicit.txt
  else
    "$MM" create --yes --no-rc --strict-channel-priority --prefix "$ENV" --file "$HERE/membrane-environment.yml"
  fi
fi
export PATH="$ENV/bin:$PATH"
"$ENV/bin/python" - <<'PY'
import openmm, parmed, MDAnalysis, shutil
from openff.toolkit import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from openmm.app import ForceField
assert shutil.which('antechamber') and shutil.which('parmchk2')
ForceField('amber14-all.xml','amber14/tip3p.xml')
print('Membrane dependencies imported. Scientific build/simulation not run.')
PY
"$MM" list --no-rc --prefix "$ENV" --explicit > /opt/yuanhe/evidence/membrane-linux-64.explicit.txt
"$ENV/bin/python" -m pip freeze > /opt/yuanhe/evidence/membrane.freeze.txt
cp "$HERE/membrane-environment.yml" /opt/yuanhe/evidence/membrane-environment.yml
echo 'Membrane environment ready. Preserve /opt/yuanhe/evidence for exact platform dependency reconstruction.'
