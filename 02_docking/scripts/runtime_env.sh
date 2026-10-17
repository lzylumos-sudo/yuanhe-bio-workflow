#!/bin/bash
export LD_LIBRARY_PATH="$(find /opt/sbvs/env/lib/python3.11/site-packages/nvidia -mindepth 2 -maxdepth 2 -type d -name lib -printf '%p:' 2>/dev/null)/opt/sbvs/env/lib:/usr/lib/wsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENMM_CPU_THREADS=2
