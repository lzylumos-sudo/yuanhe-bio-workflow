# 18-month research and commercialization plan

Planning horizon: month 0 is approval/start, not a retroactive promise. Two full-time contributors (scientific/ML and platform), part-time medicinal chemistry, access to a Linux GPU and an experimental collaborator are assumed. Actual staffing, cash budget and partner availability are unconfirmed. With one part-time developer, prioritize research reliability and one differentiation experiment; expect a longer timeline and defer enterprise services. These are targets, not commitments or revenue forecasts.

| Phase | Deliverables | Exit evidence | Owner / dependencies | Adjustment rule |
|---|---|---|---|---|
| A: months 0-3 | Rights inventory; approved branding; stable install; MurJ acceptance/export; benchmark runner; approved read-only demo | Current release gates resolved or explicitly failed; fresh reproduction by another person; actual environment/weight rights recorded | Platform + scientific lead; recovered artifacts or approved rerun; rights-holder input | Do not label stable/general-purpose until reproducible. Unclear model rights block commercial delivery |
| B: months 3-6 | 3-5-target pilot; performance/recovery fixes; multi-objective report; usability study; 8-12 targeted customer/research interviews | Held-out control results, failure rates/costs, paired user-task evidence; documented recurring customer problem | Scientific lead + chemist; external researchers and licensed labels | If no repeatable problem/value signal, narrow segment before adding infrastructure |
| C: months 6-12 | Structure-constraint ablation; fixed-pool AL; private enterprise prototype only after approval; 1-2 scoped research pilots | Independent fixed-budget baseline comparison; calibrated uncertainty/OOD failures; tenant isolation tests before customer data | ML + platform; rights-cleared data, approved private workspace, pilot partners | No method benefit: retain integration product, revise or stop algorithm direction. No customer need: defer API/queue expansion |
| D: months 12-18 | Private deployment package; support/security/backup plan; commercial terms reviewed; prospective assay collaboration; paid-pilot discussions | Signed scope and actual customer acceptance; measured support/compute margin; prospective assay outcomes if available | Founder + counsel + engineering + experimental partner | No paid intent or prohibitive support costs: revise delivery/service model; no enterprise launch by calendar alone |

## Initial product and business hypothesis

Research Edition supplies source-available, reproducible workflows for permitted research uses. The first commercial offering, if rights and customer demand are cleared, could be a scoped private deployment/service with auditable candidate-prioritization evidence, rather than a general drug-discovery subscription. Enterprise licensing, access control and support are not implemented services today. Do not promise a drug, assay success, regulatory compliance or guaranteed savings.

Interview computational groups, biotech project leads and CRO teams separately. Test who pays, what data can leave their environment, existing tools, turnaround requirements, integration burden and willingness to pay for a defined deliverable. Record negative interviews. A letter of interest is not revenue, and a paid pilot is not repeatable product-market fit. No contact messages or customer commitments are authorized by this document.

## Resource and economics model

Before Phase A execution, estimate each proposed run from a short approved pilot: cost = GPU hours x actual quoted rate + CPU/storage/egress + human review/support + failed/repeated runs + external assay/chemistry + license/legal costs. Record currency/date and separate marginal run cost from platform development. Freeze caps for model calls, GPU hours, storage and assay spend; stop on cap, do not silently upscale. No prices or expenses are assumed here.

Measure gross contribution for each pilot after direct compute, support, third-party and experimental costs. Track installation success, job failure/recovery, turnaround, human review minutes, evidence completeness, repeat use and actual payment. Stars/downloads are discovery signals, not scientific or commercial validation.

## IP and release governance

Before publishing a new method, dataset or candidate series, record inventors/owners, institution obligations, public disclosures and intended jurisdictions; obtain qualified advice and an explicit disclosure decision. Keep future proprietary implementation and customer data private under approved access. Do not defer patent review until month 12 if a disclosure is imminent. Keep existing history and grants intact. Every distribution needs a version-specific rights and SBOM review.

## Immediate order

First close release evidence and license unknowns. Next build the benchmark/acceptance runner and multi-objective evidence table. Then run the fixed-pool active-learning sandbox. Start original large-model work only if smaller baselines and data rights justify its cost. See [three priorities and feasibility](technical_strategy.md).
