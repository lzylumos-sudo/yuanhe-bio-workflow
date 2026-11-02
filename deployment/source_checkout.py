"""Refuse lossy vendor checkouts on case-insensitive filesystems."""
from pathlib import Path
import subprocess
import tempfile


def case_sensitive_directory(destination):
    with tempfile.TemporaryDirectory(prefix='.case-probe-', dir=destination) as scratch:
        probe = Path(scratch) / 'YuanheCaseProbe'
        probe.write_bytes(b'case-sensitive checkout probe')
        return not probe.with_name('yuanhecaseprobe').exists()


def checked_checkout(destination, revision):
    destination = Path(destination)
    raw = subprocess.check_output(['git', '-C', str(destination), 'ls-tree', '-r',
                                   '--name-only', '-z', revision])
    paths = raw.decode('utf-8').rstrip('\0').split('\0') if raw else []
    seen, collisions = {}, []
    for name in paths:
        # Include directory components, not only colliding complete filenames.
        for index in range(1, len(Path(name).parts) + 1):
            part = '/'.join(Path(name).parts[:index]); key = part.casefold()
            if key in seen and seen[key] != part:
                collisions.append((seen[key], part))
            seen[key] = part
    if collisions:
        if not case_sensitive_directory(destination):
            raise ValueError('Vendor contains case-colliding paths: ' + repr(collisions[:4]) +
                             '. Use a WSL ext4 directory, or enable NTFS case sensitivity on the '
                             'new destination before restoring. No checkout was attempted.')
    subprocess.run(['git', '-C', str(destination), 'checkout', '--detach', revision], check=True, timeout=900)
    # Git status alone can miss case-folded duplicates; check every tracked blob.
    entries = subprocess.check_output(['git', '-C', str(destination), 'ls-tree', '-r', '-z', revision])
    for entry in entries.split(b'\0'):
        if not entry:
            continue
        metadata, name = entry.split(b'\t', 1)
        mode, kind, expected = metadata.split()
        if kind != b'blob':
            continue
        path = destination / name.decode('utf-8')
        if mode == b'120000':
            import os
            data = os.readlink(path).encode('utf-8')
        else:
            data = path.read_bytes()
        import hashlib
        actual = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if actual != expected.decode():
            raise ValueError('Restored source blob mismatch (case sensitivity or checkout filter): ' + str(path))
