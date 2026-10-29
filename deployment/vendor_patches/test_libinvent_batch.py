"""Real CPU inference regression against the pinned official LibInvent prior."""
import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'01_molecule_generation'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--executable', required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    from src.adapters.reinvent4 import run
    import torch
    from reinvent.models.libinvent.models.decorator import Encoder
    torch.set_num_threads(1)
    for layers in [1, 2]:
        for batch in [1, 2]:
            encoder = Encoder(layers, 8, 16, 0)
            output, (h, c) = encoder(torch.tensor([[2,3,1]]*batch), torch.tensor([3]*batch))
            assert output.shape == (batch,3,8) and h.shape == c.shape == (layers,batch,8)
    rows = []
    config = json.loads((ROOT/'deployment/model_dependencies.json').read_text())['reinvent4']
    gen = ROOT/'01_molecule_generation'
    for budget in [1, 2]:
        result = run(scaffold_smiles='c1ccccc1[*:1]', checkpoint=gen/config['prior'],
            checkpoint_sha256=config['prior_sha256'], output_dir=args.out/('batch_'+str(budget)),
            budget=budget, seed=17, executable=args.executable, repo=gen/'vendor/REINVENT4', device='cpu',
            environment={'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
        assert result['report']['completed'] == budget
        manifest = json.loads((args.out/('batch_'+str(budget))/'manifest.json').read_text())
        assert manifest['commands'][-1]['returncode'] == 0
        rows.append({'requested':budget, 'completed':result['report']['completed'], 'report':result['report']})
    # This fixed public fixture produces one chemically invalid sample. The
    # inference must still finish, while the strict adapter must reject it.
    rejected_out = args.out/'invalid_batch_1'
    try:
        run(scaffold_smiles='c1ccccc1[*:1]', checkpoint=gen/config['prior'],
            checkpoint_sha256=config['prior_sha256'], output_dir=rejected_out,
            budget=1, seed=41, executable=args.executable, repo=gen/'vendor/REINVENT4', device='cpu',
            environment={'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    except RuntimeError as error:
        assert 'no valid unique molecules retaining the scaffold core' in str(error), str(error)
    else:
        raise AssertionError('The known invalid sample was accepted')
    rejected_manifest = json.loads((rejected_out/'manifest.json').read_text())
    assert rejected_manifest['status'] == 'failed'
    assert rejected_manifest['commands'][-1]['returncode'] == 0
    with (rejected_out/'raw_sampling.csv').open() as handle:
        assert len(list(csv.DictReader(handle))) == 1
    assert not (rejected_out/'candidates.json').exists()
    proof={'status':'passed','invalid_sample_rejected':True,'scope':'real official-prior CPU inference; not biological validation',
           'torch_version':torch.__version__, 'batches':rows, 'checkpoint_sha256':config['prior_sha256']}
    (args.out/'acceptance.json').write_text(json.dumps(proof,indent=2))
    print(json.dumps(proof,indent=2))


if __name__ == '__main__':
    main()

