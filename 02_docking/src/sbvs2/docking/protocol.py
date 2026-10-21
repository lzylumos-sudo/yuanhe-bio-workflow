"""Docking protocols (plan §17). Presets are [P]-class starting points, frozen into a protocol hash."""
from __future__ import annotations
from dataclasses import dataclass, asdict
from ..hashing import hash_obj

PRESETS = {
    "quick":    dict(engine="vina",  exhaustiveness=8,  num_modes=9, seeds=(11,), cnn_scoring="none"),
    "standard": dict(engine="vina", exhaustiveness=32, num_modes=9, seeds=(11,), cnn_scoring="none"),
    "topk":     dict(engine="vina", exhaustiveness=32, num_modes=9, seeds=(11, 23, 47), cnn_scoring="none"),
}
LEGACY_MODEL = "crossdock_default2018_1_3"   # legacy protocol kept; new models never mixed into old rankings


@dataclass(frozen=True)
class Protocol:
    name: str
    engine: str
    exhaustiveness: int
    num_modes: int
    seeds: tuple
    cnn_scoring: str
    cnn_model: str | None = None
    allow_cpu_fallback: bool = False
    engine_version: str = "1.2.7"
    scoring: str = "vina"
    energy_range: float = 3.0

    @classmethod
    def preset(cls, name: str, **kw):
        return cls(name=name, **{**PRESETS[name], **kw})

    def hash(self):
        return hash_obj(asdict(self))

    def gnina_args(self, receptor, ligand, box, seed, out_sdf, log):
        a = ["--receptor", receptor, "--ligand", ligand, "--center_x", box.center[0], "--center_y", box.center[1],
             "--center_z", box.center[2], "--size_x", box.size[0], "--size_y", box.size[1], "--size_z", box.size[2],
             "--exhaustiveness", self.exhaustiveness, "--num_modes", self.num_modes, "--seed", seed,
             "--cnn", self.cnn_model or LEGACY_MODEL, "--cnn_scoring", self.cnn_scoring, "--out", out_sdf, "--log", log]
        return [str(x) for x in a]    # argv array, never a shell string

    def vina_args(self, receptor_pdbqt, ligand_pdbqt, box, seed, out_pdbqt):
        a = ["--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt, "--center_x", box.center[0], "--center_y", box.center[1],
             "--center_z", box.center[2], "--size_x", box.size[0], "--size_y", box.size[1], "--size_z", box.size[2],
             "--exhaustiveness", self.exhaustiveness, "--num_modes", self.num_modes, "--seed", seed, "--scoring", self.scoring, "--energy_range", self.energy_range, "--out", out_pdbqt]
        return [str(x) for x in a]
