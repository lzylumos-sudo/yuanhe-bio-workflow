# Historical MurJ trial: receptor ready, docking blocked

Target: EcMurJ / 9NU4 chain A / UniProt P0AF16. The historical murj2 receptor had 3,785 atoms, 498 residues and author range 3-500. Removed fusion residues -101..2 comprised 104 residues and 813 atoms; component deletions were recorded. The original complex with its 29 modeled Sgl residues remained archived. Numbering round trips, selected-domain coverage and reference retention passed; WT restoration was not evaluated.

At that time no reviewed local docking box was available. 9NU4 contains a peptide reference, not a co-crystal small molecule; the upstream generic 1DRF box was inapplicable. The system therefore did not fabricate a blind box or MurJ poses. Approximately 20 chemically complete input compounds were available, but no target docking submission was made. GNINA was also blocked by missing CUDA runtime libraries.

Pose QA, local mechanism, redocking and membrane pose assessment remained not_evaluated; candidate count was zero. The required next inputs were a sourced/reviewed MurJ site hypothesis and an executable GNINA runtime. The later workbench can calculate and preview a candidate box from selected residues, but that still requires scientific review and is not retrospective validation of this historical trial.
