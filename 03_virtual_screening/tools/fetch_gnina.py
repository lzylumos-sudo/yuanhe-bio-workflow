"""Fetch the exact GNINA build recorded by the upstream project's vendor lock.

The upstream MemPilot project pins gnina 1.3.3 with a sha256 in
``env/vendor.lock.json``.  This script downloads that same artifact and refuses
to keep it unless the digest matches, so the screening studio never runs an
unverified docking binary.
"""

from __future__ import annotations

import hashlib
import json
import os
import ssl
import sys
import urllib.request

URL = "https://github.com/gnina/gnina/releases/download/v1.3.3/gnina.cuda12.8.static"
SHA256 = "3340c1f49cd3c7c84d8699182a1c6af13c7fa2a22448d1204640446106f72172"
TARGET = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vendor", "gnina")


def digest(path: str) -> str:
    sha = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            sha.update(block)
    return sha.hexdigest()


def main() -> int:
    os.makedirs(os.path.dirname(TARGET), exist_ok=True)
    if os.path.exists(TARGET) and digest(TARGET) == SHA256:
        print(json.dumps({"status": "already_present", "path": TARGET, "sha256": SHA256}))
        return 0
    context = ssl.create_default_context()
    request = urllib.request.Request(URL, headers={"User-Agent": "zvss/0.9"})
    partial = TARGET + ".part"
    with urllib.request.urlopen(request, timeout=120, context=context) as response, open(partial, "wb") as handle:
        total = 0
        while True:
            block = response.read(1 << 20)
            if not block:
                break
            handle.write(block)
            total += len(block)
    actual = digest(partial)
    if actual != SHA256:
        os.remove(partial)
        print(json.dumps({"status": "digest_mismatch", "expected": SHA256, "actual": actual}))
        return 2
    os.replace(partial, TARGET)
    os.chmod(TARGET, 0o755)
    print(json.dumps({"status": "downloaded", "path": TARGET, "bytes": total, "sha256": actual}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
