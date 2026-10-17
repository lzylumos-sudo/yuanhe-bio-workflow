"""List cache records, so stale 'unsupported' evidence can be recomputed."""

from __future__ import annotations

import json
import os
import sqlite3

PROJECT = r"C:\Users\user\Desktop\Yuanhe Bio\06_Virtual_Screening_Studio"
DB = os.path.join(PROJECT, "runtime", "zvss.sqlite3")

con = sqlite3.connect(DB)
con.row_factory = sqlite3.Row
print("total cache rows:", con.execute("SELECT COUNT(*) AS n FROM cache").fetchone()["n"])
rows = con.execute(
    "SELECT stage, completed, COUNT(*) AS n FROM cache GROUP BY stage, completed ORDER BY stage"
).fetchall()
for row in rows:
    print(" ", dict(row))
print("--- rows mentioning posebusters/unsupported in meta or path ---")
for row in con.execute("SELECT cache_key, stage, scope_id, output_path, meta FROM cache"):
    text = json.dumps({k: row[k] for k in row.keys()}, default=str)
    if "unsupported" in text.lower() or "posebusters" in text.lower():
        print(" ", dict(row))
