#!/bin/bash
# Main-text Table 1 (time and peak memory): minimap2, certified mode and heuristic mode in one series.
# Same executable (tag revision-2026), 32 threads, index built from the FASTA inside the timed run, SAM to
# /dev/null, three replicates with the three configurations interleaved; the table reports the medians.
# Then checks that the certified mode's SAM equals minimap2's (header @PG line excluded) on every data set.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"
D=results/table1_modes; mkdir -p $D
T=$D/timing.tsv; I=$D/identity.tsv
echo -e "dataset\tconfig\trep\twall_s\tmaxrss_kb" > $T
echo -e "dataset\tcertified_sam_equals_minimap2" > $I
DS=("chr22_hifi map-hifi chr22_named.fa hifi_chr22.fq" "chr22_ont map-ont chr22_named.fa ont_r10_chr22.fq"
    "chr22_clr map-pb chr22_named.fa clr_chr22.fq" "chr14_hifi map-hifi ref_chr14.fa hifi_chr14.fq"
    "chr14_ont map-ont ref_chr14.fa ont_r10_chr14.fq" "chr14_clr map-pb ref_chr14.fa clr_chr14.fq")
env_of(){ case $1 in minimap2) echo "MM2_GEO=0";; certified) echo "MM2_GEO=1";; heuristic) echo "MM2_GEO=1 MM2_GEO_FAST=1";; esac; }
for s in "${DS[@]}"; do set -- $s
  for rep in 1 2 3; do for c in minimap2 certified heuristic; do
    /usr/bin/time -f "%e %M" -o $D/t env $(env_of $c) "$MM2GEO" -ax $2 -t32 $3 $4 > /dev/null 2>/dev/null
    echo -e "$1\t$c\t$rep\t$(cat $D/t | tr ' ' '\t')" >> $T
  done; done
done
for s in "${DS[@]}"; do set -- $s
  a=$(MM2_GEO=0 "$MM2GEO" -ax $2 -t32 $3 $4 2>/dev/null | grep -v '^@PG' | md5sum | cut -c1-32)
  b=$(MM2_GEO=1 "$MM2GEO" -ax $2 -t32 $3 $4 2>/dev/null | grep -v '^@PG' | md5sum | cut -c1-32)
  echo -e "$1\t$([ "$a" = "$b" ] && echo yes || echo NO)" >> $I
done
awk -F'\t' 'NR>1{k=$1"\t"$2; t[k]=t[k]" "$4; r[k]=r[k]" "$5}
  function med(s,  a,n){n=split(s,a," "); asort(a); return a[int((n+1)/2)]}
  END{for(k in t) printf "%s\t%.0f s\t%.1f GB\n", k, med(t[k]), med(r[k])/1048576}' $T | sort > $D/summary.txt
cat $I >> $D/summary.txt; echo DONE >> $D/summary.txt
