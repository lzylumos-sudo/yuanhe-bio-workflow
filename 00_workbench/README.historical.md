# Historical WEDGE-Design scaffold method

This English reference describes the earlier rule-based prototype, not the current Peptide2Mol/REINVENT4 backend. Numbers below are historical recorded outputs, not newly executed experiments. Geometric scores and peptide overlap are hypotheses, not proof of inhibition, permeability, resistance avoidance or synthetic feasibility.

The source structures are 9NU4, 9NU5 and 9NU8. The recorded literature reference is Li et al., *Convergent MurJ flippase inhibition by phage lysis proteins*, Nature 652, 1274-1280 (2026), DOI 10.1038/s41586-026-10163-w. Consult the original publication before using its claims; this translation does not independently verify it.

## Scope and inputs

Three phage lysis peptides motivated a common-interface hypothesis: occupy the TM2/TM7 groove and central cavity to interfere with MurJ conformational closure. The prototype represents scaffolds as rigid unit chains connecting pharmacophore anchors. Substituent design and true docking are separate downstream tasks.

| Structure | Peptide | Modeled peptide range | Resolution |
|---|---|---|---|
| 9NU4 | SglM | 2-30, 29 residues | 3.6 A |
| 9NU5 | SglPP7 | 27-55, 29 residues | 3.7 A |
| 9NU8 | SglCJ3 | 64-103, 40 residues | 3.6 A |

The historical parser uses MurJ residues 5-500 (496 residues), excluding 106/7/14 extra construct residues, including fusion segments and the KYLL linker at 1-4. This historical range differs from the current configurable workbench default; inspect the actual construct rather than treating either as a new biological validation.

| Measurement | Historical value |
|---|---|
| Receptor interface residues | 22 / 36 / 40 |
| Peptide interface residues | 18 / 19 / 28 |
| Peptide contact heavy atoms | 64 / 111 / 121 |
| Contact fraction of modeled peptide | 62% / 66% / 70% |
| Scaffold heavy atoms | 20-27, median 24 |
| Common-contact residues | 14, geometric convergence |
| Literature conserved residues | 13, a distinct sequence-conservation claim |

Sequence conservation and common geometric contacts are not interchangeable. The historical evidence layer checked 24 priors: 22 PASS, one FLAG (SglM D18 to K46, 6.67 A) and one WEAK. Measurements are recalculated by `pepwedge/evidence.py`; quotations and executable check rules live in `priors/murj_priors.json`; `verify.py` reports disagreements explicitly. Terminal carboxylate geometry must consider both O and OXT, not OXT alone.

## Modules and artifacts

| Module | Operation | Historical artifact |
|---|---|---|
| M0 | Coordinate measurements and prior checks | Evidence report and seven data files |
| M1 | Contacts, burial and prior-weighted hotspots | 52 interface residues |
| M2 | Cavity voxels and consensus pharmacophores | 5,440 A^3, four subpockets, 12 features |
| M3 | Length-constrained unit enumeration | 20,307 scaffolds |
| M4 | Rigid placement and contact scoring | 3,378 scored scaffolds |
| M5 | Ranking and stratified redundancy reduction | Top 200 table and diverse SDF |

```bash
python run_evidence.py --out evidence
python run_pipeline.py --out results
python run_filters.py --out results
python tests/test_smoke.py
```

Historical timings were about 30 seconds for evidence, 50 seconds for design and 81 seconds for filtering on one CPU. These are not performance guarantees. Old SDF chemistry must pass independent RDKit readback before reuse; the v2 implementation rejects invalid historical chemical inputs.

## Interface and hotspot calculation

Contacts use heavy-atom separation below 4.5 A. Salt-bridge geometry uses Lys NZ / Arg NE,NH1,NH2 and acidic or terminal carboxyl oxygens below 4.0 A. Shrake-Rupley SASA uses 92 sphere points and a 1.4 A probe; burial is free-minus-bound SASA.

`score_i = 0.30*c_i + 0.25*dSASA_i/max(dSASA) + 0.15*n_i/max(n) + 0.20*e_i + 0.10*r_i`

Here `c` is the fraction of three peptides contacting the residue, `n` counts atom pairs, and `e/r` encode functional/resistance priors. Prior-weighted rankings do not independently validate those priors.

| Residue | Score | Convergence | dSASA, A^2 | Prior |
|---|---|---|---|---|
| L50 | 0.736 | 1.00 | 106.9 | L50W resistance |
| Q244 | 0.702 | 1.00 | 70.0 | Q244P resistance |
| I248 | 0.685 | 1.00 | 78.8 | None |
| T251 | 0.655 | 1.00 | 84.1 | None |
| I237 | 0.650 | 1.00 | 63.0 | None |
| K46 | 0.621 | 0.67 | 64.9 | Functional |
| R270 | 0.589 | 0.33 | 70.9 | Functional |
| R53 | 0.588 | 0.67 | 62.0 | E273 association |
| T113 | 0.587 | 1.00 | 52.6 | None |
| I252 | 0.587 | 1.00 | 64.3 | None |

The 14 common contacts are L50, I54, V112, T113, G234, I237, L238, V240, S241, Q244, I245, I248, T251 and I252. Positive anchor priors are K46, R53, R270 and K368. These lists define hypotheses, not mandatory interactions for every valid inhibitor.

## Pocket and pharmacophores

Voxel spacing is 0.8 A, with a peptide envelope padded by 7 A. Keep points at least 3.0 A from receptor heavy atoms. Burial uses 50 approximately uniform rays, a 1.8 A probe and 13 A cutoff; retain at least 50% occluded directions. Assign subpockets by nearest lining-residue centroid.

| Subpocket | Volume, A^3 | Mean burial | Lining residues |
|---|---|---|---|
| S1 basic bowl | 2426 | 0.73 | K46 R270 D269 Q244 S243 Y267 Y266 |
| S3 TM2/TM7 groove | 1591 | 0.60 | F42 L47 L50 L51 I54 V112 T113 P115 F117 |
| S4 TM7 floor | 720 | 0.61 | I237 L238 V240 S241 I245 L247 I248 T251 I252 S255 |
| S2 carboxylate gate | 704 | 0.94 | R53 K368 E273 V367 E57 L364 |

Features derive from acidic oxygen midpoints, aromatic centroids, polar functional atoms, terminal O/OXT midpoint and hydrophobic side-chain centroids. Align shared receptor C-alpha atoms by Kabsch to the SglCJ3 frame. Cluster same-type features within 3.0 A, with one contribution per structure. Weight: `w = w0*(0.6+0.4*c)`.

| Feature | Type | Coordinates, A | Convergence | Weight | Source |
|---|---|---|---|---|---|
| A1 | Anionic | 162.2, 166.2, 185.5 | 0.67 | 1.13 | PP7 V55 / CJ3 L103 termini |
| R2 | Aromatic | 162.6, 180.0, 188.7 | 0.67 | 1.08 | PP7 Y36 / CJ3 Y76 |
| A7 | Anionic | 174.2, 168.5, 190.6 | 0.33 | 0.81 | CJ3 E96 |
| R4 | Aromatic | 165.1, 159.2, 189.0 | 0.33 | 0.95 | PP7 W53 |
| D6 | Polar | 164.9, 179.2, 192.3 | 0.67 | 0.82 | CJ3 T80 / PP7 S40 |

| Pair | Historical distance, A |
|---|---|
| A1-A7 | 13.3 |
| A1-R2 | 14.1 |
| A7-R2 | 16.5 |
| K46-R270 | approximately 9 |
| R53-K368 | approximately 12.6 |

## Scaffold enumeration

Scaffolds follow `head-link-ring-link-...-tail`; axial span is the sum of unit spans. Precompute sorted linker-span sums and use interval lookup for `abs(total_span-target_distance) <= 0.55 A`. This reduces lookup cost, while enumerating all returned combinations still costs proportional to their number.

| Unit | Span, A |
|---|---|
| 1,4-phenylene | 2.80 |
| 1,4-cyclohexylene | 2.80 |
| Thiophene-2,5-diyl | 2.26 |
| Direct bond | 1.48 |
| Ethynyl | 4.06 |
| Methylene | 2.96 |
| Ether oxygen | 2.72 |
| Amide | 4.25 |
| Ethylene | 4.49 |
| COOH / CH2COOH head | 0 / 1.52 |
| Phenyl, pyridyl, pyrimidinyl, cyclohexyl tail | 1.40 |
| tert-Butyl / benzoic-acid tail | 0 / 4.28 |

Historical filters constrain the first linker, limit repeated linker types to two and amides to one, and require a rigid ring core. Coordinates and SMILES come from the same unit definitions. These construction rules do not establish synthesis feasibility.

```python
for feature1, feature2 in design_axes:
    for head, tail in compatible_caps:
        for rings in ring_combinations:
            need = target_span - head.span - tail.span - sum(r.span for r in rings)
            for links in interval_lookup(need, tol=0.55):
                if legal(head, links, tail):
                    scaffold = compose(head, rings, links, tail)
```

## Placement and scoring

Rotate the local +x axis onto the feature-pair axis using Rodrigues rotation, sample 12 axial angles and two mirrors, and translate the head feature to the first anchor. These 24 rigid placements are not flexible docking.

| Score term | Geometric rule | Weight |
|---|---|---|
| Clash | Squared overlap below 0.78*(r_i+r_j) | -6.0 |
| Salt bridge | Anionic carboxyl oxygen to positive anchor below 4.2 A | +3.0 |
| Hydrogen-bond proxy | N/O pairs, 2.2-3.65 A | +0.5 |
| Hydrophobic | C/S pairs below 4.6 A | +0.1 |
| Exposure | No receptor atom within 5.5 A | -0.4 |
| Buried fraction | One minus exposed fraction | +0.8 |

Ether and carbonyl oxygens do not receive anionic salt-bridge credit. Peptide overlap is the fraction of scaffold atoms within 2.3 A of any peptide atom. It measures spatial overlap, not mechanistic proof.

`J = E_contact - 1.5*P_prop + 0.6*R + 0.4*(w_F1+w_F2) + 1.6*O_peptide`

`R = min(1, D/[3*(1+RotB)])`. Rank by J, retain the top 200 and limit related head/ring-count/tail combinations for diversity. SDF poses require chemical readback before downstream use.

| Funnel criterion | Historical count |
|---|---|
| Enumerated | 20307 |
| Property penalty <= 0.75 | 3378 |
| Clash < 0.25 | 1282 |
| Peptide overlap >= 0.5 | 2924 |
| Clash < 0.25, property penalty <= 0.35, overlap >= 0.5 | 538 |

| Structural filter | Rule | Passing among 538 |
|---|---|---|
| F1 | Contact K46/R53/R270/K368 with anionic oxygen | 538 |
| F2 | No acidic oxygen separation below 3.6 A | 538 |
| F3 | Bridge polar S1/S2 and hydrophobic S3/S4 | 462 |
| F4 | Re-place in each structure; best-score range <= 5 | 32 |

The top historical scaffold failed F4: scores 9.97, 8.45 and 0.76, with a missing salt bridge and clash 0.84 in SglM. Reusing a pose without re-search can change a score from +10 to -42 because side-chain geometry differs. Retain this negative result; it motivates ensemble validation rather than deletion.

## Examples and handoff

| Example | SMILES | MW | cLogP | TPSA |
|---|---|---|---|---|
| Biphenyl pyridyl amide | `OC(=O)c1ccc(c2ccc(C(=O)Nc3cnccc3)cc2)cc1` | 318.3 | 3.2 | 79.3 |
| Biphenyl pyrimidyl amide | `OC(=O)c1ccc(c2ccc(C(=O)Nc3nccnc3)cc2)cc1` | 319.3 | 2.6 | 92.2 |
| Thiophene analogue | `OC(=O)c1ccc(c2sc(C(=O)Nc3ccccc3)cc2)cc1` | 323.4 | 4.3 | 94.6 |
| tert-Butyl amide | `OC(=O)c1ccc(c2ccc(C(=O)NC(C)(C)C)cc2)cc1` | 297.4 | 4.6 | 66.4 |
| Dianionic ether | `OC(=O)c1ccc(Oc2ccc(Oc3ccc(C(=O)O)cc3)cc2)cc1` | 350.3 | 3.9 | 93.0 |

| Artifact | Purpose |
|---|---|
| `scaffolds_ranked.csv` | Chemistry, properties, scores and contact residues |
| `diverse_scaffolds.sdf` | Coordinate-bearing hypotheses requiring chemistry validation |
| `pharmacophore.json` | Features, tolerances and weights |
| `pocket.json` | Cavity and subpocket geometry |
| `interface_hotspots.csv` | Prior-weighted interface analysis |

Historical proposals included scaffold decoration, pocket-conditioned generation, ensemble docking, state-selectivity comparisons, 100-500 ns MD and free-energy follow-up. These proposals were not completed by this prototype. The current backend uses Peptide2Mol pocket-conditioned seeds and reviewed REINVENT4 cores; it does not claim peptide-pharmacophore preservation or unimplemented inpainting.

## Parameters, limitations and validation proposals

| Parameter | Historical setting |
|---|---|
| Contact / salt cutoff | 4.5 / 4.0 A |
| SASA probe / points | 1.4 A / 92 |
| Hotspot weights | 0.30 / 0.25 / 0.15 / 0.20 / 0.10 |
| Voxel / receptor clearance / envelope padding | 0.8 / 3.0 / 7.0 A |
| Burial rays / probe / cutoff / threshold | 50 / 1.8 A / 13 A / 0.50 |
| Feature merge radius | 3.0 A |
| Length tolerance | 0.55 A |
| Enumerated ring count | 2-3 |
| Enumeration mass limit | CLI default 520 Da; pipeline passes 500 Da |
| Property envelope | MW 260-470; cLogP 0.8-4.5; TPSA <=95; HBD <=2; HBA <=9; RotB <=6; rings <=4 |
| Property penalties | MW excess /40, lower shortfall /150; cLogP upper *0.6, lower *0.4; TPSA above 90 /12 |

The historical property model uses fragment estimates, not a validated permeability model. The score has no explicit membrane, water, entropy or reliable free energy. Resolution limits side-chain certainty; efflux, metabolism, synthesis and cellular access require independent assessment.

Proposed controls include reference-peptide geometry checks, randomized feature-distance controls, mutation sensitivity, outward/inward state comparisons and prospective biochemical/cellular assays. No proposed control, MIC target, candidate mechanism, synthesis route or assay availability is presented here as experimentally established. Current release acceptance is defined in the root validation guide, separately from this historical method.

- Historical values require matching source coordinates.
- Proposed controls remain unexecuted unless separately documented.
