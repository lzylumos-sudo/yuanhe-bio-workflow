# CLI

From the repository root (no new installation system):

```bash
python yuanhe version
python yuanhe doctor
python yuanhe modules
python yuanhe protocols
python yuanhe validate --config task.json
python yuanhe run --config task.json
python yuanhe status YH-YYYYMMDD-UUID
```

`./yuanhe` also works with an available Python. `PYTHONPATH=00_workbench python -m yuanhe` is equivalent. `run` creates and validates a task only; it does not automatically launch scientific stages. Continue with the existing reviewed UI or `python 00_workbench/run_research.py start --run-id ID --stage prepare`. A protocol override is shown as custom; use a new task for frozen parameter changes.

`doctor` checks registry mappings and environment prerequisites without launching molecular inference, docking or GPU simulation. It is not a scientific validation run. A missing model/dependency remains unavailable.

Existing `run_research.py`, `python -m sbvs.cli`, `python -m sbvs2`, `python -m zvss`, ZLIS and MD entrypoints remain supported as documented in their components. New docking imports use `from docking import run_vina_job, Protocol, Box`. Model access can use `src.adapters.registry.generate(name, **kwargs)`; Peptide2Mol and REINVENT4 remain distinct algorithms.
