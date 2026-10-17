import os
import tempfile

targets = [
    os.environ["TEMP"],
    os.path.join(os.environ["TEMP"], "zvss-scratch"),
    os.path.join(os.environ["LOCALAPPDATA"], "Temp"),
]
for sub in sorted(os.listdir(os.environ["TEMP"]))[:400]:
    path = os.path.join(os.environ["TEMP"], sub)
    if os.path.isdir(path) and sub.startswith("dsh"):
        targets.append(path)
        targets.append(os.path.join(path, "zvss-scratch"))
for cand in targets:
    try:
        os.makedirs(cand, exist_ok=True)
        probe = os.path.join(cand, ".probe")
        with open(probe, "w", encoding="utf-8") as handle:
            handle.write("ok")
        os.remove(probe)
        print("WRITABLE", cand)
    except Exception as exc:  # noqa: BLE001
        print("denied  ", cand, type(exc).__name__, exc)
print("gettempdir:", tempfile.gettempdir())
