"""Extract CUDA runtime libraries from NVIDIA's pip wheels for the GNINA build.

GNINA 1.3.3 is linked against the CUDA 12 runtime; the WSL2 driver provides
``libcuda`` but not ``libcudart``.  The pip wheels published by NVIDIA are the
smallest auditable source of those shared objects, so this script resolves the
current wheel URL from the PyPI JSON API, records its sha256, and unpacks the
``lib*`` members into ``vendor/<package>/lib``.
"""

from __future__ import annotations

import hashlib
import json
import os
import ssl
import sys
import urllib.request
import zipfile

VENDOR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vendor")
PACKAGES = ["nvidia-cuda-runtime-cu12", "nvidia-cublas-cu12"]


def http_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "zvss/0.9"})
    with urllib.request.urlopen(request, timeout=120, context=ssl.create_default_context()) as response:
        return json.loads(response.read().decode("utf-8"))


def http_bytes(url: str, destination: str) -> str:
    digest = hashlib.sha256()
    request = urllib.request.Request(url, headers={"User-Agent": "zvss/0.9"})
    with urllib.request.urlopen(request, timeout=300, context=ssl.create_default_context()) as response, open(
        destination, "wb"
    ) as handle:
        while True:
            block = response.read(1 << 20)
            if not block:
                break
            digest.update(block)
            handle.write(block)
    return digest.hexdigest()


def pick_wheel(items) -> dict | None:
    if isinstance(items, dict):
        items = items.get("urls", [])
    for item in items or []:
        name = item.get("filename", "")
        if name.endswith(".whl") and "manylinux" in name and "x86_64" in name:
            return item
    return None


def main() -> int:
    record: dict[str, object] = {"packages": {}}
    for package in PACKAGES:
        meta = http_json(f"https://pypi.org/pypi/{package}/json")
        version = meta["info"]["version"]
        wheel = pick_wheel(meta["releases"].get(version) or [])
        if not wheel:
            candidates = [v for v in meta["releases"] if meta["releases"][v]]
            for candidate in sorted(candidates, reverse=True):
                wheel = pick_wheel(meta["releases"][candidate])
                if wheel:
                    version = candidate
                    break
        if not wheel:
            record["packages"][package] = {"status": "no_linux_wheel"}
            continue
        filename = wheel["filename"]
        target_dir = os.path.join(VENDOR, package, "lib")
        os.makedirs(target_dir, exist_ok=True)
        archive = os.path.join(VENDOR, filename)
        if not (os.path.exists(archive) and os.path.getsize(archive) == wheel.get("size")):
            sha = http_bytes(wheel["url"], archive)
        else:
            sha = hashlib.sha256(open(archive, "rb").read()).hexdigest()
        extracted = []
        with zipfile.ZipFile(archive) as bundle:
            for member in bundle.namelist():
                base = os.path.basename(member)
                if base.startswith("lib") and ".so" in base:
                    with bundle.open(member) as source, open(os.path.join(target_dir, base), "wb") as sink:
                        sink.write(source.read())
                    extracted.append(base)
        record["packages"][package] = {
            "version": version,
            "wheel": filename,
            "sha256": sha,
            "extracted": sorted(extracted),
        }
    with open(os.path.join(VENDOR, "cuda_runtime.json"), "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(json.dumps(record, indent=1)[:2000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
