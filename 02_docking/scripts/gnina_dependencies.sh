#!/bin/bash
set -euo pipefail
export PIP_CACHE_DIR=/opt/sbvs/cache/pip
/opt/sbvs/env/bin/python -m pip install nvidia-cudnn-cu12==9.10.2.21 nvidia-cublas-cu12==12.9.2.10 nvidia-cuda-nvrtc-cu12==12.9.86 nvidia-cuda-runtime-cu12==12.8.90 nvidia-cusparse-cu12==12.5.8.93 nvidia-cufft-cu12==11.3.3.83 nvidia-cusolver-cu12==11.7.3.90 nvidia-nvjitlink-cu12==12.8.93
