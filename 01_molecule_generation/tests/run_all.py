#!/usr/bin/env python3
"""Run all v2 unit/contract tests in the chemistry environment."""
from __future__ import annotations
import importlib.util
import pathlib
import sys
import traceback

ROOT = pathlib.Path(__file__).resolve().parent
PROJECT = ROOT.parent
sys.path.insert(0, str(PROJECT))


def main() -> int:
    failures = []
    files = sorted(ROOT.glob("test_*.py"))
    if not files:
        print("no v2 tests found", file=sys.stderr)
        return 2
    for path in files:
        print(f"\n[{path.name}]")
        spec = importlib.util.spec_from_file_location(path.stem, path)
        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)
            if hasattr(module, "run_tests"):
                module.run_tests()
            elif hasattr(module, "main"):
                code = module.main()
                if code not in (None, 0):
                    raise RuntimeError(f"returned {code}")
            else:
                import unittest
                suite = unittest.defaultTestLoader.loadTestsFromModule(module)
                result = unittest.TextTestRunner(verbosity=2).run(suite)
                if not result.wasSuccessful():
                    failures.append(path.name)
        except Exception:
            traceback.print_exc()
            failures.append(path.name)
    if failures:
        print(f"FAILED ({len(failures)}): {failures}")
        return 1
    print(f"\nall {len(files)} v2 test modules passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
