#!/bin/bash
# Supplementary Table S19, full-coverage HiFi row: stock minimap2 (MM2_GEO unset) on every primary read of the GIAB
# HG002 HiFi alignment (the union of the 25 partitions of script 24), counted in streaming without keeping the output.
# Mapped reads: reads with >= 1 tp:A:P record; non-secondary records: tp:A:P records (primary and supplementary).
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/hifi53x_alignment_rate; mkdir -p $D
IDX=grch38_hifi.mmi; [ -s $IDX ] || MM2_GEO=0 "$MM2GEO" -x map-hifi -d $IDX "$GRCH38_REF" > /dev/null 2>&1
[ -s "$HG002_HIFI_BAM" ] || { echo "run: bash scripts/fetch_data.sh genomewide"; exit 1; }
"$SAMTOOLS" view -@8 -u -F 0x900 "$HG002_HIFI_BAM" | "$SAMTOOLS" fastq -@8 - 2>/dev/null \
 | tee >(awk 'NR%4==2{n++; b+=length($0)} END{printf "input reads %d, bases %d\n", n, b}' > $D/input.txt) \
 | MM2_GEO=0 "$MM2GEO" -cx map-hifi -t"$THREADS" $IDX - 2>/dev/null \
 | awk -F'\t' '/tp:A:P/{p++; if($1!=last){m++; last=$1}} {r++} END{printf "PAF records %d, non-secondary records %d, mapped reads %d\n", r, p, m}' > $D/paf.txt
sleep 5; cat $D/input.txt $D/paf.txt | tee $D/summary.txt
