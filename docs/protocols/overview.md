# Scientific protocols

The registry distinguishes `murj-discovery@1.0.0`, `soluble-protein-md@1.0.0` and `murj-membrane-md@1.0.0`. These are new registrations of audited source recipes from commit e4906f2; they do not assert that historical runs used these IDs. Historical quality labels such as `single-500ns-quality/2` remain unchanged and are not product versions.

The registry freezes source hashes, research configuration and the soluble YAML / membrane MDP recipe. `check_governance.py` detects recipe/default drift. `check_scientific_baseline.py` protects 127 scientific/config/data/history files byte-for-byte. Changes to force fields, electrostatics, state preparation, integration, thresholds or selection semantics require an explicit new protocol revision and validation scope. UI-only changes do not.

See [MurJ discovery](murj_discovery.md), [soluble MD](soluble_md.md), [membrane MD](membrane_md.md), and existing [method and dependency notes](../15_VINA_RESP_RESEARCH.md). The registry's engineering recipe is not a universally literature-recommended protocol and does not establish efficacy, convergence or equilibrium.
