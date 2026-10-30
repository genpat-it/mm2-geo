#!/bin/bash
# E. coli ONT timing (Supplementary section on the E. coli spike-in): reads SRR9900640 against the spike-in
# reference, -ax map-ont, 32 threads, SAM to /dev/null; minimap2 (MM2_GEO=0), certified mode (MM2_GEO=1) and
# heuristic mode (MM2_GEO=1 MM2_GEO_FAST=1), same executable (tag revision-2026), one warm-up run, then three
# interleaved replicates. Output: results/ecoli_timing/{timing.tsv,summary.txt} (wall s, peak RSS in kB / GiB).
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/ecoli_timing; mkdir -p $D
REF=ecoli_vc/ref_mut.fa; RD=ecoli_real2/SRR9900640.fastq
E_mm2="MM2_GEO=0"; E_cert="MM2_GEO=1"; E_heur="MM2_GEO=1 MM2_GEO_FAST=1"
MM2_GEO=0 "$MM2GEO" -ax map-ont -t32 $REF $RD > /dev/null 2>&1   # warm-up (page cache)
echo -e "config\trep\twall_s\tmaxrss_kB" > $D/timing.tsv
for rep in 1 2 3; do for c in mm2 cert heur; do v=E_$c
  /usr/bin/time -f "%e %M" -o $D/t env ${!v} "$MM2GEO" -ax map-ont -t32 $REF $RD > /dev/null 2>/dev/null
  echo -e "$c\t$rep\t$(tr ' ' '\t' < $D/t)" >> $D/timing.tsv; done; done; rm -f $D/t
awk -F'\t' 'NR>1{t[$1]=t[$1]" "$3; r[$1]=r[$1]" "$4}
  function med(s, a,n){n=split(s,a," "); asort(a); return a[int((n+1)/2)]}
  END{for(k in t) printf "%s\t%.1f s\t%.1f GiB\n", k, med(t[k]), med(r[k])/1048576}' $D/timing.tsv | sort > $D/summary.txt
cat $D/summary.txt
