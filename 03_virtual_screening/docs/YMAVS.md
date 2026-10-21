# Y-MAVS — Yuanhe Membrane-Aware Virtual Screening Engine

Status: implementation complete for the P0 membrane modules and the P1 ensemble data
layer. This document states exactly what is implemented, what was executed, and what
remains a scientific or data limitation.

## What Y-MAVS changes

Y-MAVS is a post-GNINA evidence layer inside `zvss`. It does not re-implement docking
and does not replace GNINA. It adds:

1. membrane region classification from the existing verified PPM3 orientation transfer;
2. ligand-defined membrane pocket analysis (depth, residue composition, annotated TM
   segment contact);
3. static geometric accessibility from a probe-expanded receptor grid;
4. membrane-compatibility descriptors and a van der Waals overlap proxy;
5. optional per-conformational-state evidence integration without averaging scores;
6. an evidence-first ranking that preserves the original GNINA values and baseline rank;
7. additive qualification requirements for membrane geometry and pocket localization.

The system is **not** a trained binding-affinity model. `membrane_aware_score` is an
ordinal priority inside the current candidate pool, not ΔG or activity probability.

## Modules and call graph

```text
stages.run_all
  └─ membrane.environment.run
       ├─ membrane.orientation.transfer          (existing, reused)
       ├─ membrane.regions.MembraneGeometry      (new)
       ├─ membrane.pocket.MembranePocketAnalyzer (new)
       ├─ membrane.accessibility.MembraneAccessibilityAnalyzer (new)
       ├─ membrane.compatibility.MembraneCompatibilityAnalyzer (new)
       ├─ membrane.ensemble.MembraneEnsembleAnalyzer (optional; new)
       └─ outputs under runs/<run_id>/stages/membrane/
  └─ scoring.aggregate.run
       └─ scoring.membrane_aware.MembraneAwareScorer (new)
  └─ select.pipeline.run
       ├─ select.qualify (extended, additive only)
       └─ scoring.membrane_aware.write_ranking_files
            └─ reports/membrane_aware_ranking.csv/json,
               membrane_evidence.json, membrane_screening_report.html
  └─ md.package (unchanged; uses the selected pose_id)
```

`membrane` still runs before `scoring`, so qualification reads measured evidence rather
than a placeholder. A target opts in with:

```yaml
screening_mode: membrane_aware
membrane_screening:
  enabled: true
```

A filename or a PDB code is never used to infer membrane mode. `generic.yaml` remains
`screening_mode: generic` and follows the old path exactly.

## Coordinate and evidence rules

* The existing `orientation.transfer` is reused. A PPM3 plane is accepted only after
  receptor SHA, atom coverage and rigid-fit RMSD checks. No new PCA guess is made when
  sourced evidence exists.
* If no orientation evidence exists, a PCA plane from the docking receptor may be used
  only as `proxy_only`; its thickness is the configured `geometry.proxy_thickness_A`,
  not a measured bilayer. It cannot satisfy `membrane_geometry` formally.
* Every pose must belong to the same receptor SHA as the plane. A mismatch is recorded
  as `needs_input` and the membrane sub-analyses are not used for ranking or formal
  qualification.
* Positive/negative depth means only a side of the supplied normal. Y-MAVS never calls
  a side cytoplasmic or extracellular without independent topology evidence.
* Membrane core is a depth range, not a lipid-occupied volume or a whole bilayer.

## Qualification semantics

The existing statuses are retained:

`pass`, `fail`, `not_evaluated`, `needs_input`, `unsupported`.

`select.qualify` adds three optional requirement ids:

* `membrane_geometry` — verified plane, finite pose depth, receptor-hash binding;
* `pocket_localization` — ligand-defined pocket contacts and a known membrane region;
* `accessibility_evidence` and `ensemble_evidence` — optional evidence requirements.

The target block can add requirements but cannot remove the budget's base requirements:

```yaml
membrane_screening:
  qualification:
    formal_requirements: [membrane_geometry, pocket_localization]
    accessibility_evidence: auxiliary   # or required
    ensemble_evidence: optional         # or required
```

Accessibility is auxiliary by default. A static no-path result is `uncertain`, not a
failure. Formal qualification remains separate from ranking.

## Outputs

Written under the run directory, never the repository root:

```text
runs/<run_id>/stages/membrane/
  membrane_geometry.json
  residue_membrane_regions.csv
  ligand_membrane_regions.csv
  pose_membrane_evidence.json
  pocket_analysis.json
  pocket_residues.csv
  accessibility_analysis.json
  accessibility_analysis.csv
  accessibility_paths.json
  membrane_compatibility.json
  membrane_compatibility.csv
  ensemble_analysis.json
  ensemble_state_coverage.csv

runs/<run_id>/reports/
  membrane_aware_ranking.csv
  membrane_aware_ranking.json
  membrane_evidence.json
  membrane_screening_report.html
```

`membrane_aware_ranking.csv` preserves `compound_id`, `state_id`, `pose_id`,
`baseline_rank`, `membrane_aware_rank`, `rank_shift`, raw GNINA affinity/CNN values,
membrane descriptors, evidence status, formal qualification and ranking reasons. The
HTML report is a backend-generated review artifact, not a frontend rewrite.

## Ensemble evidence

Multiple conformational states are read from per-state JSON manifests with real scores.
A state that has no manifest remains `missing`; it is never assigned zero. The schema is:

```json
{
  "schema": "ymavs-ensemble-state-evidence/1",
  "structure_id": "9NU5",
  "receptor_sha256": "<sha256 of the state receptor>",
  "coordinate_frame": "state_docking_receptor",
  "compounds": {
    "CMP...": {
      "state_id": "STA...",
      "pose_id": "POS...",
      "gnina_affinity": -8.4,
      "cnn_score": 0.71,
      "cnn_affinity": 7.2,
      "pose_qc_status": "pass",
      "state_specific_contacts": 18,
      "state_specific_accessibility": {"aqueous": "accessible", "membrane_lateral": "uncertain"}
    }
  }
}
```

Each configured state also needs its own structure/hash and its own membrane orientation
evidence. The primary state's orientation is never copied to another state. A C-alpha
Kabsch fit checks chain, residue identity, missing residues and RMSD; it does not claim
that conformational states are experimentally interchangeable.

## Usage

Install the module's declared dependencies in the existing chemistry environment, then:

```bash
zvss --root . run \
  --input /path/to/gnina/results \
  --target configs/targets/ecmurj.yaml \
  --budget configs/budgets/full.yaml
```

For another membrane protein, copy `configs/targets/membrane_protein_template.yaml` and
provide:

* a docking receptor and its hash;
* a sourced OPM/PPM3 placement or explicitly defensible same-frame orientation;
* a real GNINA result set imported by the existing path.

No docking, MD, membrane construction or ligand generation is launched by this module.

## Validation performed

Lightweight checks only:

* synthetic membrane depth/region assertions;
* synthetic pocket and accessibility grids, including a real connected path;
* focused scoring order test proving the baseline rank and raw GNINA values survive;
* proxy-orientation qualification test proving it cannot satisfy `membrane_geometry`;
* ensemble test proving missing states are not averaged or zero-filled, plus a C-alpha
  Kabsch alignment test on a rotated/translated synthetic state;
* an offline full-pipeline run using the archived 1DRF GNINA fixture with a deliberately
  proxy-only membrane target. The run completed, produced all membrane artifacts, kept
  the original GNINA scores and pose ID, and produced **zero** formal candidates because
  the membrane geometry was proxy-only.
* the real PPM3 transfer and region classification were exercised against the 9NU4
  evidence and its prepared receptor. A membrane-only run on one archived 9NU4 pose
  produced a verified plane, 13 pocket contacts, `TM2`/`TM7` annotation contact,
  membrane-embedded pocket classification and an evaluated accessibility result.
  The archived full-complex and prepared receptor share identical coordinates for all
  498 matched C-alpha atoms; that fixture normalization is recorded in the validation
  summary rather than hidden.

Not claimed: MD execution, a new docking campaign, retention of the PPM3 frame through a
different receptor preparation, benchmark enrichment (EF1%, BEDROC), or improved hit rate.
Pose reselection is deliberately disabled (`allow_pose_reselection: false`) because moving
the representative pose would invalidate pose-bound strain, redocking and qualification
evidence; per-pose membrane evidence is retained for a future version that recomputes
those checks for the new pose. Pocket volume is `not_evaluated` unless a validated volume
method is added.

## Next scientific validation

1. Run Y-MAVS on a receptor preparation that exactly matches the PPM3 evidence and retain
   the SHA chain in the report.
2. Compare plain GNINA ordering and Y-MAVS ordering on a labelled active/inactive set using
   EF1%, BEDROC and confidence intervals, not a handful of anecdotes.
3. Establish whether membrane localization/accessibility strata are stable across
   conformational states and receptor preparations.
4. Prospectively test a prespecified subset before claiming improved hit finding.
5. Only then consider a calibrated or learned ranking weight; the current layer is
   uncalibrated and intentionally rule-based.
