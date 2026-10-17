"""Quick composition probe for the 9NU4 receptor question (analysis only)."""

from __future__ import annotations

import collections
import json
import os
import sys

from zvss import compat

compat.bootstrap()

from zvss.util import parse_pdb_text, read_text  # noqa: E402

ref = sys.argv[1]
name = sys.argv[2] if len(sys.argv) > 2 else "structure"
st = parse_pdb_text(read_text(os.path.join(ref, name)), keep_hydrogens=False)
prot = [a for a in st.atoms if not a.record.startswith("HETATM")]
het = [a for a in st.atoms if a.record.startswith("HETATM")]
print("==", name, "total", len(st.atoms), "protein", len(prot), "hetero", len(het))
print("  chains", dict(collections.Counter(a.chain for a in prot)))
print("  hetero", dict(collections.Counter((a.resname, a.chain) for a in het)))

by_chain = collections.defaultdict(lambda: collections.defaultdict(int))
for a in prot:
    by_chain[a.chain][a.resseq] += 1
for chain, residues in sorted(by_chain.items()):
    nums = sorted(residues)
    print(f"  chain {chain}: {len(nums)} residues, auth {nums[0]}..{nums[-1]}, atoms {sum(residues.values())}")
    below3 = [n for n in nums if n < 3]
    if below3:
        print(f"    auth<3 residues: {len(below3)} ({below3[0]}..{below3[-1]}), atoms {sum(residues[n] for n in below3)}")
    over500 = [n for n in nums if n > 500]
    if over500:
        print(f"    auth>500 residues: {len(over500)} ({over500[0]}..{over500[-1]}), atoms {sum(residues[n] for n in over500)}")
    native = [n for n in nums if 3 <= n <= 500]
    print(f"    native 3..500 residues: {len(native)}, atoms {sum(residues[n] for n in native)}")

mapping_path = os.path.join(ref, "mapping_9NU4_A.json")
if os.path.exists(mapping_path):
    with open(mapping_path, encoding="utf-8") as handle:
        mapping = json.load(handle)
    rows = mapping.get("rows") or []
    print("  mapping rows", len(rows), "keys", list(mapping))
    mapped = [r for r in rows if r.get("uniprot_resseq")]
    print("  mapped", len(mapped), "first", mapped[0] if mapped else None)
