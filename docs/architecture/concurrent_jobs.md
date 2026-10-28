# Concurrent generation and receptor preparation

The workbench permits one pipeline stage and one CPU receptor-preparation
operation at the same time, in the same run or in different runs. Complex
preparation must already be complete before preparing a receptor. The receptor
operation writes to its own immutable operation directory. Applying its reviewed
result changes the future docking review; it never changes a running generator's
frozen input or overwrites a receptor already used by docking.

## Resource and lifecycle contract

| Work | Capacity | Isolation |
|---|---|---|
| Pipeline stage, including generation | One per installation | Existing GPU allocation retained; one pipeline writer per run |
| Receptor preparation | One per installation | CPU only, at most two library threads, unique operation directory |
| Export, environment operations, MD tools, diagnostics, source/evidence revisions | Exclusive | Wait until background work exits |
| Settings for an unstarted stage | Independent | Per-stage lock; dispatched/completed inputs remain immutable |

An occupied slot returns a specific capacity error. There is no automatic queue
or automatic retry. This conservative two-job policy does not promise arbitrary
multi-GPU throughput, CPU/memory reservations or a cluster scheduler. Existing
GPU UUID leases still protect the scientific workers.

Every job has its own process group, identity (PID and Linux start time), status,
log and output. **Stop this operation** affects only the selected operation;
**Stop this stage** affects the selected pipeline stage; **Stop all jobs in this
run** affects all active jobs in that run and leaves other runs alone. Termination
allows checkpoint flushing before escalation. A cancelled record does not free
its execution slot until the worker and inheriting descendants release it.

All lock descriptors are inherited by engine subprocesses. The original global
lock is retained as a shared execution/exclusive maintenance barrier, including
compatibility with an older server holding its exclusive lock. A server restart
does not grant a second writer access to a live worker's files. Status checks
match individual stage/operation records, so an active generator cannot conceal
an abandoned receptor operation.

## Reference and acceptance boundary

This design follows the publicly documented principles of asynchronous jobs,
resource limits and dependency checks in Schrödinger
[Jobcontrol](https://learn.schrodinger.com/public/python_api/2025-3/jobcontrol.html)
and [JobDJ](https://learn.schrodinger.com/public/python_api/2026-1/api/schrodinger.job.queue.html).
It does not use Schrödinger implementation code and does not assert Suite
compatibility, feature parity, certification or equivalent scientific accuracy.

`00_workbench/tests/test_concurrent_jobs.py` tests real Linux processes, both
launch orders, shared/cross-run execution, capacity limits, inherited locks,
failure cleanup, service reattachment, input immutability and cancellation.
The controlled workers in these tests do not perform model inference.

`deployment/validate_concurrent_receptor.py` additionally runs the actual
PDBFixer/OpenMM receptor worker on a public PDB while a controlled long-lived
pipeline process holds the generation slot. It verifies overlapping lifetimes,
CPU isolation, preparation output, preserved input hashes and successful review
application. This is an orchestration/receptor engineering check, not generation,
docking, affinity, MD convergence or scientific validation.

Before wider commercial deployment, separately accept durable queues and restart
reconciliation, administrator resource budgets, multi-user authorization and data
isolation, operational backups/monitoring, dependency/model usage rights, and
target-specific scientific protocols. These remain separate release gates.

## Upgrade

Deploy the merged source through the administrator workflow. Do not restart
running scientific workers to apply this patch. Frozen historical runs retain
their source checks and may require a separately reviewed source revision or a
new run; this patch does not silently relabel their provenance. GitHub merging
alone does not update an already running workstation or a previous release ZIP.
