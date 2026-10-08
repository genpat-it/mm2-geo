#!/bin/bash
# Band-doubling rounds and joint (d,E) statistics per gap (Supplementary Table S15).
# Instrumented build (MM2_GEO_GAPLOG); heuristic mode and stock, chr22 HiFi/ONT/CLR and the genome-wide subset.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; mkdir -p results/gaplog
H="$(cd "$(dirname "$0")" && pwd)"
for s in "hifi map-hifi chr22_named.fa hifi_chr22.fq" "ont map-ont chr22_named.fa ont_r10_chr22.fq" "clr map-pb chr22_named.fa clr_chr22.fq" "wg138k map-hifi $GRCH38_REF $GW_READS"; do set -- $s
  for g in 1 0; do
    MM2_GEO=$g MM2_GEO_GAPLOG=results/gaplog/${1}_geo$g.tsv "$MM2GEO_INSTR" -cx $2 -t"$THREADS" $3 $4 > results/gaplog/${1}_geo$g.paf 2>/dev/null
    gzip -f results/gaplog/${1}_geo$g.tsv
  done
done
python3 "$H/18_gap_stats.py" results/gaplog | tee results/gap_rounds_dE.txt
