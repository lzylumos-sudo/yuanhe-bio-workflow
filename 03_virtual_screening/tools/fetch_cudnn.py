"""Fetch a cuDNN 9 runtime archive from NVIDIA's public redist CDN.

GNINA 1.3.3 (cuda12.8 static build) dynamically links ``libcudnn.so.9``.  WSL2
ships the GPU driver but not cuDNN, so the runtime has to be provided next to
the binary.  The archive is taken from NVIDIA's official redistributable CDN and
its sha256 is recorded, matching how the GNINA binary itself is pinned.
"""

from __future__ import annotations

import hashlib
import json
import os
import ssl
import sys
import tarfile
import urllib.request

BASE = "https://developer.download.nvidia.com/compute/cudnn/redist/cudnn/linux-x86_64/"
PREFERRED = "cudnn-linux-x86_64-9.8.0.87_cuda12-archive.tar.xz"
VENDOR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vendor")


def digest(path: str) -> str:
    sha = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            sha.update(block)
    return sha.hexdigest()


def list_dir() -> list[str]:
    context = ssl.create_default_context()
    request = urllib.request.Request(BASE, headers={"User-Agent": "zvss/0.9"})
    with urllib.request.urlopen(request, timeout=120, context=context) as response:
        html = response.read().decode("utf-8", "replace")
    names = []
    for part in html.split("href=")[1:]:
        quote = part[0]
        end = part.find(quote, 1)
        if end > 0:
            names.append(part[1:end])
    return sorted(name for name in names if name.endswith(".tar.xz") and "_cuda12" in name)


def main() -> int:
    names = list_dir()
    print("available cuda12 archives:", names[-5:])
    name = PREFERRED if PREFERRED in names else names[-1]
    url = BASE + name
    archive = os.path.join(VENDOR, name)
    os.makedirs(VENDOR, exist_ok=True)
    if os.path.exists(archive) and os.path.getsize(archive) > 1 << 20:
        print("archive already present:", archive)
    else:
        print("downloading", url)
        context = ssl.create_default_context()
        request = urllib.request.Request(url, headers={"User-Agent": "zvss/0.9"})
        with urllib.request.urlopen(request, timeout=300, context=context) as response, open(
            archive + ".part", "wb"
        ) as handle:
            while True:
                block = response.read(1 << 20)
                if not block:
                    break
                handle.write(block)
        os.replace(archive + ".part", archive)
    sha = digest(archive)
    target = os.path.join(VENDOR, "cudnn")
    os.makedirs(target, exist_ok=True)
    with tarfile.open(archive, "r:xz") as tar:
        members = [m for m in tar.getmembers() if "/lib/" in m.name or "/include/" in m.name]
        for member in members:
            member.name = member.name.split("/", 1)[1]
        tar.extractall(target, members=members)
    libs = sorted(os.listdir(os.path.join(target, "lib")))
    record = {
        "archive": os.path.basename(archive),
        "sha256": sha,
        "url": url,
        "libcudnn_present": "libcudnn.so.9" in libs,
        "libs_sample": libs[:10],
    }
    with open(os.path.join(VENDOR, "cudnn_runtime.json"), "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(json.dumps(record, indent=1))
    return 0 if record["libcudnn_present"] else 2


if __name__ == "__main__":
    sys.exit(main())
