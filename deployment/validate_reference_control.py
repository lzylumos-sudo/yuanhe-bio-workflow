"""Real fixed public Vina reference acceptance; no candidate qualification is emitted."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'00_workbench'))
from scientific.validation import docking_control


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vina', required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args(); out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    result = docking_control(out, {'vina':str(Path(args.vina).resolve()), 'chem_python':sys.executable}, {}, out, ROOT)
    assert len(result['runs']) == 3 and result['candidate_qualified'] is False
    assert not list(out.rglob('*.qualification.json'))
    for folder in (out/'standalone/jobs').iterdir():
        proof = json.loads((folder/'reference_control.json').read_text())
        rec = json.loads((folder/'job_record.json').read_text())
        assert proof['schema'] == 'yuanhe.reference-control/1'
        assert rec['status'] == 'ok' and not rec['candidate_qualified']
    (out/'acceptance.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'status':'completed', 'seeds':len(result['runs']),
                     'top1_recovery_fraction':result['top1_recovery_fraction'],
                     'candidate_qualified':False}))


if __name__ == '__main__':
    main()
