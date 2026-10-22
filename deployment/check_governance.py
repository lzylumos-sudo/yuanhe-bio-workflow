#!/usr/bin/env python3
"""Fail closed on registry, protocol or generated software-version drift."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'00_workbench'))
from yuanhe.registry import validate
from yuanhe.protocols import validate_sources
from yuanhe.version import identity
import json
print(json.dumps({'registry':validate(),'protocols':validate_sources(),'software':identity()},indent=2))
