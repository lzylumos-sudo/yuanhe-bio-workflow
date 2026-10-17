import copy
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from modules.utils.embedded_chemistry import EmbeddedChemistryError, PREFIX, read_embedded, write_embedded

coordinate = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C  \nEND\n"
manifest = {"system_type": "protein_ligand", "ligand": {"resname": "LIG", "chain": "L", "resid": "1", "atom_mapping": {"C1": 1}}}
# The container is chemistry-agnostic; production consumers validate the actual SDF.
sdf = "container_roundtrip_fixture\n$$$$\n"
wrapped = write_embedded(coordinate, manifest, sdf)
assert all(len(line) <= 80 for line in wrapped.splitlines() if line.startswith(PREFIX))
assert "".join(line for line in wrapped.splitlines(keepends=True) if not line.startswith(PREFIX)) == coordinate
payload = read_embedded(wrapped)
assert payload["manifest"] == manifest and payload["ligand_sdf"] == sdf
assert read_embedded(write_embedded(coordinate, {"system_type": "protein_protein"}))["ligand_sdf"] is None
assert read_embedded(coordinate) is None
assert write_embedded(wrapped, manifest, sdf) == wrapped
lines = wrapped.splitlines(keepends=True)
cases = ["".join(line for line in lines if "DATA 000001 " not in line),
         wrapped + next(line for line in lines if "DATA 000001 " in line),
         wrapped.replace("SHA256 1 ", "SHA256 3 ", 1)]
for bad in cases:
    try:
        read_embedded(bad)
    except EmbeddedChemistryError:
        pass
    else:
        raise AssertionError("Corrupt metadata was accepted")
print("EMBEDDED_CHEMISTRY_TEST_PASS: roundtrip, coordinate preservation, <=80 columns, missing/duplicate/hash corruption rejection")
