# MurJ discovery protocol

Purpose: traceable structure-informed candidate generation, qualification, Vina docking, ZVSS/Y-MAVS evidence filtering, MD and reviewed experimental handoff.

Authority: `00_workbench/configs/research_dual_blackwell.json` and registered source hash. Structure is a task input (default project 9NU4); its bytes and chain/range/pocket settings are frozen in the run. Peptide2Mol is pocket-conditioned inference. REINVENT4 uses LibInvent prior scaffold decoration; this path does not execute RL or fine-tuning.

Budgets and Vina settings are listed in the [configuration reference](../reference/configuration.md). Chemical identity, risk and purchase/synthesis evidence gates remain in their existing implementations. Sampling attempts, valid samples and unique parent/state identities remain different counts. Scores are empirical docking scores, not experimental affinity probabilities.

Outputs: raw model records, valid/failed/duplicate counts, qualified candidates, per-job poses, evidence decisions, system handoffs and reviewed reports. Tool/model versions come from actual execution records; unavailable models are unknown, not verified from configured hashes alone.

Validation: strict default/hash checks, existing identity/preflight/screening regressions and lightweight structure preparation. Full model sampling, real Vina, target-GPU membrane MD and wet-lab validation were not performed in this refactor. References, assumptions and original scientific limitations remain in `docs/15_VINA_RESP_RESEARCH.md` and component method documentation. The configured budgets are project choices, not universal scientific recommendations.
