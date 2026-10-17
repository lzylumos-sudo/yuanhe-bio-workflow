#!/bin/bash
# Resumable GNINA download with stall detection. Verifies sha256 against env/vendor.lock.json (legacy lock value).
cd "$(dirname "$0")/../vendor" || exit 1
U=https://github.com/gnina/gnina/releases/download/v1.3.3/gnina.cuda12.8.static
F=gnina-1.3.3.cuda12.8.static
EXP=3340c1f49cd3c7c84d8699182a1c6af13c7fa2a22448d1204640446106f72172
for i in $(seq 1 200); do
  curl -L --fail -C - --speed-limit 20000 --speed-time 30 --connect-timeout 20 -o $F.part "$U" && break
  echo "retry $i size=$(stat -c %s $F.part)"; sleep 3
done
GOT=$(sha256sum $F.part | cut -d' ' -f1)
echo "$GOT  $F" > $F.sha256
if [ "$GOT" = "$EXP" ]; then mv $F.part $F && chmod +x $F && echo VERIFIED_MATCHES_VENDOR_LOCK; else echo "HASH_MISMATCH got=$GOT expected=$EXP"; fi
