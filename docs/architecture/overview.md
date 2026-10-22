# Yuanhe Bio Platform architecture

Yuanhe Bio Platform is an AI-Driven Drug Discovery Platform. The module identity is independent of directory numbering and execution order.

| Module ID | Official name | Implementations |
|---|---|---|
| GEN | Molecular Generation | Peptide2Mol; REINVENT4 LibInvent prior sampling |
| DOCK | Molecular Docking | Official Vina; preserved sbvs/sbvs2 interfaces; GNINA history |
| SCREEN | Virtual Screening | ZVSS containing Y-MAVS membrane evidence |
| INTERACT | Interaction Analysis | ZLIS static interactions |
| MD | Molecular Dynamics | OpenMM soluble; OpenMM builder; ParmEd conversion; GROMACS membrane |
| WORKBENCH | Research Workbench | Persistent workflow, research controls, reviews and provenance |

Governance lives in `00_workbench/yuanhe`; it does not replace numerical scientific implementations. Registries live in root `configs`. Original numbered directories, scripts, environment locks and licenses remain in place.

See [module inventory](module_inventory.md), [call graph](dependency_graph.md), [duplication classifications](duplication_audit.md), [registry contract](module_registry.md) and [migration scope](../migration/architecture_refactor.md).

The workbench is the authority for integrated research tasks. Independent components retain their native runtime stores and schemas for compatibility; they must not be mistaken for a newly consolidated scheduler or database.
