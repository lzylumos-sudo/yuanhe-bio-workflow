"""Probe: assessment rows for one compound, and WSL/GNINA reachability."""

from __future__ import annotations

import os
import sqlite3
import subprocess
import sys

PROJECT = r"C:\Users\user\Desktop\Yuanhe Bio\06_Virtual_Screening_Studio"


def db_probe(compound: str) -> None:
    path = os.path.join(PROJECT, "runtime", "zvss.sqlite3")
    con = sqlite3.connect(path)
    con.row_factory = sqlite3.Row
    for row in con.execute("SELECT value FROM schema_info WHERE key='schema_version'"):
        print("schema_version:", row["value"])
    query = (
        "SELECT run_id, metric, value, value_text, status, reason_code FROM assessment "
        "WHERE scope_id=? AND stage='strain' ORDER BY run_id, metric"
    )
    for row in con.execute(query, (compound,)):
        print(" ", dict(row))


def run(argv: list[str], timeout: int = 60) -> None:
    try:
        proc = subprocess.run(argv, capture_output=True, timeout=timeout)
        out = (proc.stdout or b"").decode("utf-8", "replace").strip()
        err = (proc.stderr or b"").decode("utf-8", "replace").strip()
        print(f"$ {' '.join(argv)} -> rc={proc.returncode}")
        if out:
            print("  out:", out[:600].replace("\n", " | "))
        if err:
            print("  err:", err[:600].replace("\n", " | "))
    except Exception as exc:  # noqa: BLE001
        print(f"$ {' '.join(argv)} -> EXC {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "db":
        db_probe(sys.argv[2] if len(sys.argv) > 2 else "CMP75C7A4B0B41E4390")
    else:
        run(["wsl.exe", "--status"])
        run(["wsl.exe", "-l", "-q"])
        run(["wsl.exe", "-d", "Ubuntu-22.04", "-u", "root", "-e", "bash", "-lc", "which gnina && gnina --version | head -2"])
        run(["docker", "--version"])
        run(["ssh", "-V"])
