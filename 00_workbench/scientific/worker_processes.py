"""Keep detached engine children inside the worker's cancellation lifecycle."""
import os
from pathlib import Path
import signal
import time


def descendants(pid):
    result = []
    pending = [pid]
    while pending:
        parent = pending.pop()
        try:
            children = [int(x) for x in Path(f'/proc/{parent}/task/{parent}/children').read_text().split()]
        except OSError:
            continue
        for child in children:
            try:
                text = Path(f'/proc/{child}/stat').read_text()
                identity = text[text.rfind(')')+2:].split()[19]
            except (OSError, IndexError):
                continue
            result.append((child, identity)); pending.append(child)
    return result


def alive(pid, identity):
    try:
        text = Path(f'/proc/{pid}/stat').read_text()
        fields = text[text.rfind(')')+2:].split()
        return fields[19] == identity and fields[0] != 'Z'
    except (OSError, IndexError):
        return False


def stop_children(children, grace=5):
    for signum in [signal.SIGTERM, signal.SIGKILL]:
        for pid, identity in reversed(children):
            if alive(pid, identity):
                try:
                    os.kill(pid, signum)
                except ProcessLookupError:
                    pass
        deadline = time.monotonic()+grace
        while time.monotonic() < deadline and any(alive(*row) for row in children):
            time.sleep(.05)


def install():
    def stop(signum, frame):
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        stop_children(descendants(os.getpid()))
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM, stop)
