#!/bin/bash
# Speed of the heuristic mode under different thread counts, output handling, index scope and builds
# (Supplementary Table S20). Genome-wide HiFi subset ($GW_READS) against a chr1-only index and the full GRCh38;
# submission-2026 executable; bioconda minimap2 and an -O3 rebuild for the build comparison. Median of 3 runs.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/speedchecks; mkdir -p $D
[ -s $D/chr1.mmi ] || { "$SAMTOOLS" faidx "$GRCH38_REF" CM000663.2 > $D/chr1.fa; MM2_GEO=0 "$MM2GEO_SUB" -x map-hifi -d $D/chr1.mmi $D/chr1.fa > /dev/null 2>&1; }
[ -s grch38_hifi.mmi ] || MM2_GEO=0 "$MM2GEO" -x map-hifi -d grch38_hifi.mmi "$GRCH38_REF" > /dev/null 2>&1
head -n 120000 "$GW_READS" > $D/reads_30k.fq
run(){ # label executable env threads index reads output
  /usr/bin/time -f "%e %M" -o $D/t env $3 "$2" -ax map-hifi -t$4 $5 $6 > $7 2>/dev/null; echo -e "$1\t$(cat $D/t)" | tee -a $D/timing.tsv; rm -f $D/out.sam; }
for t in 1 4 8 32; do rd="$GW_READS"; [ $t = 1 ] && rd=$D/reads_30k.fq
  for out in devnull disk; do o=/dev/null; [ $out = disk ] && o=$D/out.sam
    for rep in 1 2 3; do for g in 0 1; do run "t$t:$out:geo=$g:rep$rep" "$MM2GEO_SUB" MM2_GEO=$g $t $D/chr1.mmi "$rd" $o; done; done; done; done
for rep in 1 2 3; do for g in 0 1; do run "grch38:t32:geo=$g:rep$rep" "$MM2GEO_SUB" MM2_GEO=$g 32 grch38_hifi.mmi "$GW_READS" /dev/null; done; done
for t in 4 32; do for rep in 1 2 3; do
  run "bioconda:t$t:rep$rep" "$BIOCONDA_MM2" MM2_GEO=0 $t $D/chr1.mmi "$GW_READS" /dev/null
  for g in 0 1; do run "O3:t$t:geo=$g:rep$rep" "$MM2GEO_O3" MM2_GEO=$g $t $D/chr1.mmi "$GW_READS" /dev/null; done; done; done
same=$(cmp -s <(MM2_GEO=1 "$MM2GEO_SUB" -cx map-hifi -t8 chr22_named.fa hifi_chr22.fq 2>/dev/null) <(MM2_GEO=1 "$MM2GEO" -cx map-hifi -t8 chr22_named.fa hifi_chr22.fq 2>/dev/null) && echo yes || echo NO)
echo "heuristic mode identical in submission-2026 and revision-2026 (chr22 HiFi): $same" | tee -a $D/summary.txt
