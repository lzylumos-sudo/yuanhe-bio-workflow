"""Find the WSL distro that actually owns /opt/sbvs/vendor/gnina."""

from __future__ import annotations

import subprocess


def decode(raw: bytes) -> str:
    for encoding in ("utf-16-le", "utf-8"):
        try:
            text = raw.decode(encoding)
        except UnicodeDecodeError:
            continue
        if "\x00" not in text:
            return text.replace("\r", "")
    return raw.decode("utf-8", "replace")


def wsl(args: list[str]) -> tuple[int, str]:
    proc = subprocess.run(["wsl.exe", *args], capture_output=True)
    text = decode(proc.stdout).strip() + ("\n[stderr] " + decode(proc.stderr).strip() if proc.stderr else "")
    return proc.returncode, text


def main() -> None:
    rc, listing = wsl(["-l", "-q"])
    names = [line.strip() for line in listing.splitlines() if line.strip() and not line.startswith("[stderr]")]
    print(f"wsl -l -q rc={rc} distros={names}")
    for name in names:
        rc, out = wsl(["-d", name, "-u", "root", "-e", "bash", "-lc",
                       "test -x /opt/sbvs/vendor/gnina && echo HAVE_GNINA; /opt/sbvs/vendor/gnina --version 2>&1 | head -2"])
        print(f"--- {name} rc={rc}\n{out}")


if __name__ == "__main__":
    main()
