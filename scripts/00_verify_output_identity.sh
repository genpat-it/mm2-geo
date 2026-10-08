#!/bin/bash
# Output identity checks (Methods): builds stock minimap2 v2.30 from https://github.com/lh3/minimap2 (tag v2.30)
# and compares, on the chr22 HiFi, ONT and CLR slices (full SAM records, @PG header line excluded):
#   A) mm2-geo with MM2_GEO unset            == stock minimap2 v2.30
#   B) certified mode (MM2_GEO=1 MM2_GEO_CERT=1) == stock minimap2 v2.30
#   C) heuristic mode, revision-2026         == heuristic mode, submission-2026
#   D) heuristic mode != stock (the feature is active)
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/identity; mkdir -p $D
STOCK="$BIN/minimap2_v2.30"
if [ ! -x "$STOCK" ]; then
  rm -rf "$BIN/src/minimap2" && git clone -q --branch v2.30 --depth 1 https://github.com/lh3/minimap2 "$BIN/src/minimap2"
  make -C "$BIN/src/minimap2" -j8 minimap2 > /dev/null && cp "$BIN/src/minimap2/minimap2" "$STOCK"
fi
h(){ env $1 "$2" -ax $3 -t"$THREADS" chr22_named.fa $4 2>/dev/null | grep -v '^@PG' | sha256sum | cut -c1-16; }
for s in "hifi map-hifi hifi_chr22.fq" "ont map-ont ont_r10_chr22.fq" "clr map-pb clr_chr22.fq"; do set -- $s
  st=$(h MM2_GEO=0 "$STOCK" $2 $3); off=$(h MM2_GEO=0 "$MM2GEO" $2 $3)
  ce=$(h "MM2_GEO=1 MM2_GEO_CERT=1" "$MM2GEO" $2 $3)
  he=$(h MM2_GEO=1 "$MM2GEO" $2 $3); hs=$(h MM2_GEO=1 "$MM2GEO_SUB" $2 $3)
  { echo "== $1"
    echo "   A) MM2_GEO unset == minimap2 v2.30:            $([ $st = $off ] && echo yes || echo NO)"
    echo "   B) certified == minimap2 v2.30:                $([ $st = $ce ] && echo yes || echo NO)"
    echo "   C) heuristic revision-2026 == submission-2026: $([ $he = $hs ] && echo yes || echo NO)"
    echo "   D) heuristic differs from minimap2 (active):   $([ $he != $st ] && echo yes || echo NO)"; } | tee -a $D/summary.txt
done
