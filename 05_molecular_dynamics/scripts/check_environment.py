#!/usr/bin/env python3
"""Diagnose the complete Ubuntu MD stack; optional tiny CPU/checkpoint smoke."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from modules.automation.backend import dependency_report
from modules.utils.core import atomic_json


def check_command(name, flag):
    executable = shutil.which(name)
    if not executable:
        return {"name": name, "passed": False, "reason": "executable missing"}
    try:
        with tempfile.TemporaryDirectory(prefix="nexus-env-check-") as folder:
            result = subprocess.run([executable, flag], capture_output=True, text=True, timeout=20, cwd=folder)
        output = (result.stdout + result.stderr).strip()
        loader_error = any(token in output.lower() for token in ("error while loading shared libraries", "cannot open shared object", "traceback (most recent call last)", "command not found", "segmentation fault"))
        return {"name": name, "passed": result.returncode in (0, 1) and bool(output) and not loader_error,
                "returncode": result.returncode, "output_tail": output[-1600:]}
    except Exception as exc:
        return {"name": name, "passed": False, "reason": f"{type(exc).__name__}: {exc}"}


def native_smoke():
    import openmm as mm
    from openmm import unit
    system = mm.System()
    system.addParticle(12)
    force = mm.CustomExternalForce("0.5*k*(x*x+y*y+z*z)")
    force.addGlobalParameter("k", 100)
    force.addParticle(0, [])
    system.addForce(force)
    integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, .002 * unit.picosecond)
    integrator.setRandomNumberSeed(20261016)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("CPU"), {"Threads": "2"})
    context.setPositions([[.1, 0, 0]])
    context.setVelocitiesToTemperature(300 * unit.kelvin, 20261016)
    integrator.step(1000)
    before = context.getState(getEnergy=True)
    energy = before.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not math.isfinite(energy):
        raise RuntimeError("Nonfinite CPU smoke-test energy.")
    checkpoint = context.createCheckpoint()
    del context, integrator
    resumed_integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, .002 * unit.picosecond)
    resumed_integrator.setRandomNumberSeed(20261016)
    resumed = mm.Context(system, resumed_integrator, mm.Platform.getPlatformByName("CPU"), {"Threads": "2"})
    resumed.loadCheckpoint(checkpoint)
    after = resumed.getState(getEnergy=True)
    loaded_ps = after.getTime().value_in_unit(unit.picosecond)
    if not math.isclose(loaded_ps, 2.0, abs_tol=1e-8):
        raise RuntimeError("Checkpoint roundtrip time differs from two ps.")
    result = {"passed": True, "platform": "CPU", "threads": 2, "steps": 1000,
              "binary_checkpoint_reloaded": True, "time_ps": float(loaded_ps), "openmm": mm.__version__,
              "scope": "One-particle installation smoke; no protein, ligand parameterization or scientific MD validation."}
    del resumed, resumed_integrator, system
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", action="store_true", help="Run a 1000-step one-particle CPU/checkpoint smoke test.")
    parser.add_argument("--output", help="Save the complete diagnostic JSON.")
    args = parser.parse_args()
    result = dependency_report()
    result["command_health"] = []
    if result["target_host"]:
        for name, flag in (("antechamber", "-h"), ("parmchk2", "-h"), ("sqm", "-h"), ("tleap", "-h"),
                           ("MMPBSA.py", "--version"), ("cpptraj", "--version"), ("sander", "-h")):
            result["command_health"].append(check_command(name, flag))
    result["ready"] = result["ready"] and all(c["passed"] for c in result["command_health"])
    if args.native and result["ready"]:
        try:
            result["native_smoke"] = native_smoke()
        except Exception as exc:
            result["native_smoke"] = {"passed": False, "error": f"{type(exc).__name__}: {exc}"}
            result["ready"] = False
    elif args.native:
        result["native_smoke"] = {"passed": False, "status": "NOT_RUN", "reason": "Target host or dependency checks failed."}
    result["status"] = "PASS" if result["ready"] and (not args.native or result["native_smoke"]["passed"]) else "NOT_READY"
    if args.output:
        atomic_json(Path(args.output), result)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
