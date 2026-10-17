# Zhiyu Ligand Interaction Studio (ZLIS) 0.1.0

Protein-ligand interaction analysis with atom-level tables, editable 2D diagrams, PyMOL scenes and real-atom PDB files. Both platform launchers use `src/zlis`.

## Start

- macOS: `Start_Mac.command` or `bash run.sh web`; CLI: `bash run.sh analyze inputs/3PTB.pdb --out results`.
- Windows: `Start_Windows.cmd` / `run.ps1`; actual Windows execution remains untested.
- Initial macOS environment: `bash env/setup_mac.sh`. Keep the Conda environment on the internal filesystem, not exFAT.
- Relocated paths with spaces and Unicode are covered by a pipeline regression.

Commands: `ligands <pdb>`; `analyze <pdb> [--ligand A:BEN1] [--template lig.sdf|SMILES] [--charge N] [--model N] [--offline] [--plip] [--no-png] [--layout layout.json]`; `web [--port 8765]`. The web service binds to 127.0.0.1 and requires its startup token; it has no account system and must not be exposed publicly.

## Outputs

Each ligand receives `interactions.csv/.xlsx/.json`, publication/detailed/audit 2D SVG/PNG/PDF, an HTML layout editor and audit view, `complex_analysis.pdb`, `ligand.sdf`, PyMOL scene/scripts, scene configuration, binding-site image and optional PLIP crosscheck. Save `layout.json` after moving labels and redraw using `--layout`.

## Scientific contracts

- ProLIF with `count=True` and explicit ZLIS geometry define interactions. PyMOL and 2D views read the resulting JSON without recomputing assignments.
- Stable atom identity is `M{model}:{chain}:{resname}{resseq}{icode}:{atomname}`; added hydrogens append `:H{n}`.
- Ligand heavy-atom coordinates always come from the input PDB. Templates supply chemistry only. Maximum heavy-atom displacement is recorded and must remain 0.0 A.
- Hydrogen-bond arrows point donor to acceptor, not electron flow. Other interactions have no arrows.
- Water bridges are candidates with both O-to-polar-atom legs within 3.5 A; roles are unspecified because crystallographic water usually lacks hydrogens. Report the longer leg. No water is invented.
- No interaction energies are calculated or fabricated.
- Default RCSB component protonation is marked `needs_confirmation`. The BEN +1 template is an explicit hypothesis.
- OpenMM adds protein hydrogens at pH 7.0; residue protonation is a model assumption. Missing OXT termini use a temporary atom that is later removed; terminal carboxyl hydrogens and charges are not modeled.

See `docs/TEST_REPORT.md` for historical checks and limitations.
