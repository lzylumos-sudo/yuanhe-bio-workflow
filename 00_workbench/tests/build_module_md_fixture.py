"""Portable real-engine MD fixture: a tiny argon system, not a biological model."""
from pathlib import Path
import argparse
import subprocess


def build(gmx, target):
    target = Path(target); target.mkdir(parents=True, exist_ok=True)
    (target/'system.top').write_text("""[ defaults ]
1 2 yes 0.5 0.8333
[ atomtypes ]
AR 18 39.948 0.0 A 0.34 0.997
[ moleculetype ]
ARG 1
[ atoms ]
1 AR 1 ARG AR 1 0.0 39.948
[ system ]
Public engineering argon fixture
[ molecules ]
ARG 2
""")
    (target/'system.gro').write_text(
        'Public engineering argon fixture\n    2\n'
        '    1ARG     AR    1   1.000   1.000   1.000\n'
        '    2ARG     AR    2   1.500   1.000   1.000\n'
        '   3.00000   3.00000   3.00000\n')
    (target/'production.mdp').write_text("""integrator = md
dt = 0.001
nsteps = 100
cutoff-scheme = Verlet
nstlist = 10
rlist = 1.0
coulombtype = Cut-off
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
tcoupl = no
pcoupl = no
gen-vel = yes
gen-temp = 20
gen-seed = 17
nstlog = 10
nstenergy = 10
nstxout-compressed = 10
comm-mode = Linear
nstcomm = 10
""")
    with (target/'grompp.log').open('w') as log:
        subprocess.run([str(gmx), 'grompp', '-f', 'production.mdp', '-c', 'system.gro',
            '-p', 'system.top', '-o', 'production.tpr'], cwd=target, stdout=log,
            stderr=subprocess.STDOUT, check=True)
    return target/'production.tpr'


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--gmx', required=True); ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args(); print(build(args.gmx, args.out))
