#!/bin/bash
# Fixed small band versus the adaptive band (Supplementary Table S8): minimap2 with -r 20,20, minimap2 with
# default settings and the heuristic mode (m=20 HiFi, m=50 ONT) on chr22; alignment time (median of 3, index
# build included, as in Table 1) and small-variant F1 (Clair3 + rtg vcfeval, GIAB v4.2.1 high-confidence regions).
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"
S=fr_scratch; mkdir -p $S; R=hifi_bench/results_FIXED_R.txt; : > $R
export LC_ALL=C
med(){ printf '%s\n' "$@" | sort -n | awk '{a[NR]=$1}END{print a[int((NR+1)/2)]}'; }
tc(){ zcat "$1" 2>/dev/null | awk '!/^#/{if(length($4)==1&&length($5)==1)s++;else i++}END{printf "%d %d",s+0,i+0}'; }
run(){ # cfg preset reads -> command
  case $1 in default) echo "env MM2_GEO=0 $MM2GEO -ax $2 -t$THREADS chr22_named.fa $3";;
             r20)     echo "env MM2_GEO=0 $MM2GEO -ax $2 -r 20,20 -t$THREADS chr22_named.fa $3";;
             geo)     echo "env MM2_GEO=1 MM2_GEO_FAST=1 MM2_GEO_MARGIN=$4 $MM2GEO -ax $2 -t$THREADS chr22_named.fa $3";; esac; }
for t in "HiFi map-hifi hifi_chr22.fq $MODEL_HIFI hifi 20" "ONT map-ont ont_r10_chr22.fq $MODEL_ONT ont 50"; do
  read tech preset reads model plat m <<< "$t"
  for cfg in default r20 geo; do
    cmd="$(run $cfg $preset $reads $m)"; ts=()
    for i in 1 2 3; do /usr/bin/time -f %e -o $S/t bash -c "$cmd > /dev/null 2>/dev/null"; ts+=("$(cat $S/t)"); done
    eval "$cmd" 2>/dev/null | "$SAMTOOLS" sort -@8 -o $S/tmp.bam - 2>/dev/null
    "$SAMTOOLS" reheader -c 'sed "s/SN:chr22/SN:CM000684.2/"' $S/tmp.bam > $S/b.bam 2>/dev/null; "$SAMTOOLS" index -@8 $S/b.bam
    conda_on "$CONDA_ENV_CLAIR3"; rm -rf $S/c3
    "$CLAIR3" --bam_fn=$S/b.bam --ref_fn=chr22_cm.fa --threads="$THREADS" --platform=$plat --model_path="$model" \
       --output=$S/c3 --ctg_name=CM000684.2 --bed_fn=conf_chr22.bed > $S/c3.log 2>&1
    conda_off; conda_on "$CONDA_ENV_RTG"; rm -rf $S/ev
    "$RTG" vcfeval -b truth_chr22.vcf.gz -c $S/c3/merge_output.vcf.gz -t ref_chr22_sdf -e conf_chr22.bed -o $S/ev > /dev/null 2>&1
    read tps tpi < <(tc $S/ev/tp.vcf.gz); read fps fpi < <(tc $S/ev/fp.vcf.gz); read fns fni < <(tc $S/ev/fn.vcf.gz); read tbs tbi < <(tc $S/ev/tp-baseline.vcf.gz)
    conda_off
    awk -v ts=$tps -v ti=$tpi -v fps=$fps -v fpi=$fpi -v fns=$fns -v fni=$fni -v tbs=$tbs -v tbi=$tbi -v l="$tech $cfg" -v t=$(med "${ts[@]}") \
     'function f(tp,fp,tb,fn, p,r){p=(tp+fp)?tp/(tp+fp):0;r=(tb+fn)?tb/(tb+fn):0;return sprintf("%.4f",(p+r)?2*p*r/(p+r):0)}
      BEGIN{printf "  %-14s time %7.1fs  SNV F1 %s  INDEL F1 %s\n",l,t,f(ts,fps,tbs,fns),f(ti,fpi,tbi,fni)}' | tee -a $R
  done
done
rm -rf $S; echo "### DONE" | tee -a $R
