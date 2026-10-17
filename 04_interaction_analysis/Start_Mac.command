#!/usr/bin/env bash
# Double-click on macOS: starts the local web UI (127.0.0.1 only).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE" || exit 1
bash "$HERE/run.sh" web
status=$?
if [ $status -ne 0 ]; then echo; echo "ZLIS exited with status $status. Press Return to close."; read -r _; fi
