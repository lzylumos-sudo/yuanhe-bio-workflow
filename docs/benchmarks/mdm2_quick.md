# MDM2–p53 retrospective benchmark

The repository CLI exposes `yuanhe benchmark`. It reuses the entity-aware receptor preparation, Meeko chemistry mapping, pinned Vina 1.2.7 adapter, risk/availability gates and existing peptide residue-contact overlap. It does not change MD defaults or other target configurations.

The research execution entry requires a frozen manifest and never writes product qualification evidence. Product Vina still fails closed on missing acquisition evidence. `current_product` runs the existing risk/acquisition rules; `ppi_benchmark` applies explicitly declared research property limits and annotates availability. An empty deliverable list is an honest outcome.

## Lifecycle

Use a Python environment containing RDKit, Meeko, molscrub, OpenMM, PDBFixer, BioPython, numpy/scipy, pandas, psutil and filelock. Set OMP/BLAS thread limits to 1; Vina has its explicit CPU budget.

```
python yuanhe benchmark preflight --root RUN_ROOT --vina /path/to/vina --started-at ISO8601_WITH_TIMEZONE
python yuanhe benchmark acquire --root RUN_ROOT
python yuanhe benchmark prepare-target --root RUN_ROOT
python yuanhe benchmark throughput --root RUN_ROOT
python yuanhe benchmark curate --root RUN_ROOT
python yuanhe benchmark prepare-library --root RUN_ROOT --count 550
python yuanhe benchmark freeze --root RUN_ROOT
python yuanhe benchmark run --root RUN_ROOT
python yuanhe benchmark evaluate --root RUN_ROOT
```

`curate` encodes the explicitly reviewed primary table mappings for this benchmark, not a generic classifier of ChEMBL records. Downloaded activity records are not automatically accepted. The provenance bundle must retain original primary sources. Receptor/ligand preparation errors and incomplete jobs remain in the frozen evaluation denominator. Only the curator and evaluator open `evaluation/`; computation reads anonymous `inputs/`.

Freeze after reviewing source/geometry and completing development throughput, before any formal test docking. Do not adjust parameters after reference recovery is observed. Successful identical jobs are reused; source/input hash changes fail closed. Existing run start/deadline is preserved by preflight. Stop starting docking 15 minutes before the 180 minute deadline and allow time for report/export. A complete run is not automatically a successful scientific validation.

A0 ranks Vina affinity. A1 ranks existing peptide residue-contact recall then affinity. A2 uses a fixed, geometrically plausible rigid perturbation of the same contact atoms. All share the identical best Vina pose; this tests post-docking peptide reranking, not peptide-conditioned search or complete chemical pharmacophore equivalence.

Top-k recall uses all frozen active parents. EF uses fixed k slots and the full frozen prevalence; unfilled slots are misses. Unknown backgrounds are not experimental negatives. 17/20 actives belong to one paper series, so broad scaffold generalization is not established. No MD, MM/PBSA or FEP runs.

Root scientific Top-50 is explicitly `ppi_benchmark/A1`. All six profile/arm lists are in `profiles/`; deliverable lists additionally require risk/acquisition qualification. Do not present research recovery as product success.

Chemical state selection preserves the standardized source tautomer and specified tetrahedral/double-bond stereochemistry, with uniform molscrub acid/base enumeration at pH 7.4. Unspecified stable stereochemistry is rejected before library freezing. All preparation failures remain in the evaluation denominator.

A chemistry-integrity failure invalidated the first batch before unblinding before any recovery metrics were opened. The amended run keeps the original start and deadline. `yuanhe benchmark reuse --root NEW_ROOT --from-root OLD_ROOT` reuses a successful job only when the actual regenerated ligand PDBQT bytes, receptor, engine, search box, protocol, seed and CPU match. Original raw job records remain unchanged; an explicit binding records changed anonymous IDs and input-file hashes. Initial development throughput is an estimate from the earlier preparation policy, not a second validation population.
