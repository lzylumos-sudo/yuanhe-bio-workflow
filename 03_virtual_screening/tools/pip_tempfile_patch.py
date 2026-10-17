"""Run pip with a patched ``tempfile``.

Why this exists
---------------
On this host the DSH file sandbox refuses writes into directories that were
created with POSIX mode ``0o700`` (``os.mkdir(path, 0o700)``).  CPython's
``tempfile`` module creates every temporary directory that way, and pip uses
``tempfile`` for its build tracker, unpack directory and ephemeral wheel
cache, so a plain ``pip install`` dies with::

    PermissionError: [Errno 13] Permission denied: '...\\pip-unpack-xxxx\\*.metadata'

This launcher patches ``tempfile`` so temporary directories are created with the
default (world-writable-under-umask) mode, then hands control to pip.  It does
not change what pip installs, only where it stages files.

Usage
-----
    python tools/pip_tempfile_patch.py --python <target-python> install ...
    python tools/pip_tempfile_patch.py install ...            # current python
"""

from __future__ import annotations

import os
import sys
import tempfile

_real_mkdir = tempfile._os.mkdir  # type: ignore[attr-defined]


def _permissive_mkdir(path, mode=0o777, *args, **kwargs):  # noqa: ANN001
    return _real_mkdir(path, 0o777, *args, **kwargs)


tempfile._os.mkdir = _permissive_mkdir  # type: ignore[attr-defined]

# ``mkdtemp`` records the first failure's errno; keep behaviour identical apart
# from the directory mode, which is all the sandbox cares about.


def main(argv: list[str]) -> int:
    from pip._internal.cli.main import main as pip_main

    return pip_main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
