#!/bin/bash
# SAM-level identity of the certified mode on the two ONT genome-wide subsets (Supplementary Tables S19, S22):
# the reads of scripts 34 (R9.4.1) and 27 (R10.4.1) aligned to the full GRCh38 with -ax map-ont, minimap2
# (MM2_GEO=0) vs certified mode (MM2_GEO=1), same executable (tag revision-2026); SAM compared without @PG.
# Scripts 34 and 27 compare PAF with CIGAR (-cx); this adds the SAM comparison (primary, secondary and
# supplementary records, SEQ/QUAL and all tags).
# Output: results/ont_genomewide_sam_identity/summary.tsv
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/ont_genomewide_sam_identity; mkdir -p $D
[ -s grch38_ont.mmi ] || MM2_GEO=0 "$MM2GEO" -x map-ont -d grch38_ont.mmi "$GRCH38_REF" > /dev/null 2>&1
echo -e "data_set\treads\tminimap2_sam_md5\tcertified_sam_md5\tcertified_sam_equals_minimap2" > $D/summary.tsv
for s in "R9.4.1 $ONT_R9_READS" "R10.4.1 ont_r10_genomewide_10pct.fq"; do set -- $s
  [ -s "$2" ] || { echo "missing $2 (run scripts 34 and 27 first)"; exit 1; }
  n=$(awk 'NR%4==1{n++}END{print n+0}' "$2")
  a=$(MM2_GEO=0 "$MM2GEO" -ax map-ont -t"$THREADS" grch38_ont.mmi "$2" 2>/dev/null | grep -v '^@PG' | md5sum | cut -c1-32)
  b=$(MM2_GEO=1 "$MM2GEO" -ax map-ont -t"$THREADS" grch38_ont.mmi "$2" 2>/dev/null | grep -v '^@PG' | md5sum | cut -c1-32)
  echo -e "$1\t$n\t$a\t$b\t$([ "$a" = "$b" ] && echo yes || echo NO)" >> $D/summary.tsv
done
cat $D/summary.tsv
