# Molecular generation and preparation

Implementation: `src/`. Standalone CLI: `murj_pipeline.py`. From this directory:

```bash
python murj_pipeline.py --run-id demo_prepare prepare
python murj_pipeline.py --run-id demo_prepare status
PYTHONPATH=. python -m pytest tests -q
```

Preparation reads `../data/structures/9nu4.pdb`. Restore models and weights into `vendor/` as described in the root setup guide. Preparation does not invoke generative models; real generation never falls back to historical rule enumeration.

`README.historical.md` records earlier development. Current entry points and acceptance scope are defined by the root documentation.
