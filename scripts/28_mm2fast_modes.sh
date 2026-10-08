#!/bin/bash
# Composition with mm2-fast (Supplementary Tables S13/S22): mm2-fast (commit 14fe36c) vs the same build with
# geo_mm2fast.patch in the heuristic and certified modes; certified identity vs mm2-fast; 3 interleaved runs.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/mm2fast; mkdir -p $D
[ -x "$MM2FAST" ] && [ -x "$MM2FASTGEO" ] || { echo "mm2-fast executables missing (see build_tools.sh)"; exit 1; }
for s in "hifi map-hifi hifi_chr22.fq" "ont map-ont ont_r10_chr22.fq" "clr map-pb clr_chr22.fq"; do set -- $s
  a=$("$MM2FAST" -cx $2 -t"$THREADS" chr22_named.fa $3 2>/dev/null | sha256sum | cut -c1-16)
  b=$(MM2_GEO=1 MM2_GEO_CERT=1 "$MM2FASTGEO" -cx $2 -t"$THREADS" chr22_named.fa $3 2>/dev/null | sha256sum | cut -c1-16)
  echo "$1: certified identical to mm2-fast: $([ $a = $b ] && echo yes || echo NO)" | tee -a $D/summary.txt
  for rep in 1 2 3; do for c in mm2fast heuristic certified; do
    case $c in mm2fast) X="$MM2FAST"; E="MM2_GEO=0";; heuristic) X="$MM2FASTGEO"; E="MM2_GEO=1";; certified) X="$MM2FASTGEO"; E="MM2_GEO=1 MM2_GEO_CERT=1";; esac
    /usr/bin/time -f "%e %M" -o $D/t env $E "$X" -ax $2 -t"$THREADS" chr22_named.fa $3 > /dev/null 2>/dev/null
    echo -e "$1\t$c\trep$rep\t$(cat $D/t)" | tee -a $D/timing.tsv; done; done
done
