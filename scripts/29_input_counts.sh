#!/bin/bash
# Input reads, bases and alignment rates per data set (Supplementary Table S19); minimap2 and heuristic mode.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results; mkdir -p $D
echo -e "dataset\tmode\tinput_reads\tinput_bases\tmapped_reads\tprimary_records\tmapped_pct" > $D/input_counts.tsv
[ -s grch38_hifi.mmi ] || MM2_GEO=0 "$MM2GEO" -x map-hifi -d grch38_hifi.mmi "$GRCH38_REF" > /dev/null 2>&1
for s in "chr22_hifi map-hifi chr22_named.fa hifi_chr22.fq" "chr22_ont map-ont chr22_named.fa ont_r10_chr22.fq" "chr22_clr map-pb chr22_named.fa clr_chr22.fq" "chr14_hifi map-hifi ref_chr14.fa hifi_chr14.fq" "chr14_ont map-ont ref_chr14.fa ont_r10_chr14.fq" "chr14_clr map-pb ref_chr14.fa clr_chr14.fq" "genomewide_subset_hifi map-hifi grch38_hifi.mmi $GW_READS"; do set -- $s
  read nr nb < <(awk 'NR%4==2{n++;b+=length($0)}END{print n, b}' $4)
  for g in 0 1; do MM2_GEO=$g "$MM2GEO" -cx $2 -t"$THREADS" $3 $4 2>/dev/null | awk -v n=$nr -v nb=$nb -v ds=$1 -v g=$g '$0~/tp:A:P/{p++; s[$1]=1} END{m=length(s); printf "%s\t%s\t%d\t%d\t%d\t%d\t%.2f\n", ds, (g?"heuristic":"minimap2"), n, nb, m, p, 100*m/n}' | tee -a $D/input_counts.tsv; done
done
