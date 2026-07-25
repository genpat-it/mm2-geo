#!/bin/bash
# Thread-scaling: alignment wall-clock vs threads, HiFi chr22, 100k-read subset. The reference
# index is built on the fly from FASTA and included in the timed stage; it is identical for all
# configs, so it cancels in the speedup ratios below.
# single run at 1-2 threads (slow), median-of-3 at >=4 threads.
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"; MF="$MM2FAST"; MFG="$MM2FASTGEO"
export LC_ALL=C
R=hifi_bench/results_THREADSCALING.txt; : > $R
S=ts_scratch; mkdir -p $S
med(){ printf '%s\n' "$@"|sort -n|awk '{a[NR]=$1}END{print a[int((NR+1)/2)]}'; }
SUB=$S/sub100k.fq; head -400000 hifi_chr22.fq > $SUB   # 100k reads
echo "subset reads: $(($(wc -l < $SUB)/4))" | tee -a $R

cmd(){ case $1 in
  mm2)     echo "MM2_GEO=0 $MG -ax map-hifi -t$2 chr22_named.fa $SUB";;
  geo)     echo "MM2_GEO=1 MM2_GEO_MARGIN=20 $MG -ax map-hifi -t$2 chr22_named.fa $SUB";;
  fast)    echo "$MF -ax map-hifi -t$2 chr22_named.fa $SUB";;
  fastgeo) echo "MM2_GEO=1 MM2_GEO_MARGIN=20 $MFG -ax map-hifi -t$2 chr22_named.fa $SUB";;
 esac; }
# runs: single at <=2 threads, median-of-3 at >=4
trun(){ local c; c="$(cmd "$1" "$2")"; local n=3; [ "$2" -le 2 ] && n=1; local t=()
  for i in $(seq 1 $n); do /usr/bin/time -f "%e" bash -c "$c >/dev/null 2>/dev/null" 2>$S/t.t; t+=("$(cat $S/t.t)"); done
  med "${t[@]}"; }

echo "### THREAD SCALING — alignment time (s, index build included), HiFi chr22 100k reads (single@1-2t, median-3@>=4t)" | tee -a $R
# speedup columns are mm2_time / cfg_time (i.e. how many times faster than mm2); >1 means faster
printf "%-4s | %-9s %-9s %-9s %-9s | %-9s %-9s\n" "thr" "mm2" "geo" "fast" "fastgeo" "mm2/geo" "mm2/fg" | tee -a $R
for thr in 1 2 4 8 16 32 64; do
  a=$(trun mm2 $thr); b=$(trun geo $thr); c=$(trun fast $thr); d=$(trun fastgeo $thr)
  awk -v t=$thr -v a=$a -v b=$b -v c=$c -v d=$d 'BEGIN{
    printf "%-4s | %-9s %-9s %-9s %-9s | %-9.2f %-9.2f\n", t, a, b, c, d, (b>0?a/b:0), (d>0?a/d:0)}' | tee -a $R
done
rm -rf $S; echo "### DONE-TS" | tee -a $R
