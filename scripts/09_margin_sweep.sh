#!/bin/bash
# Margin ablation: chr22 HiFi + ONT, m in {5,10,20,40,80}.
# time (median-3) + CIGAR-identity vs stock (all m) + clair3 F1 (m=5,20,80).
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"
SAM="$SAMTOOLS"
C3="$CLAIR3"
M_HIFI="$MODEL_HIFI"
M_ONT="$MODEL_ONT"
AR=chr22_named.fa; CR=chr22_cm.fa; CT=CM000684.2; TR=truth_chr22.vcf.gz; BD=conf_chr22.bed; SF=ref_chr22_sdf
source "$CONDA_PROFILE"
export LC_ALL=C; q(){ "$@" 2>/dev/null; }
D=marg; mkdir -p $D; R=hifi_bench/results_MARGIN.txt; : > $R
med(){ printf '%s\n' "$@"|sort -n|awk '{a[NR]=$1}END{print a[int((NR+1)/2)]}'; }
tc(){ zcat "$1" 2>/dev/null | awk '!/^#/{if(length($4)==1&&length($5)==1)s++;else i++}END{printf "%d %d",s+0,i+0}'; }
cigs(){ $SAM view -F 0x900 - 2>/dev/null | awk 'BEGIN{OFS="\t"}{print $1,$6}' | sort -k1,1; }

echo "### MARGIN ABLATION (chr22). time=median-3 32t; CIGARid vs stock; F1 clair3" | tee -a $R
mrun(){ local tech=$1 preset=$2 reads=$3 mdl=$4
  echo ">> $tech" | tee -a $R
  # stock CIGAR reference
  MM2_GEO=0 $MG -ax $preset -t32 $AR $reads 2>/dev/null | cigs > $D/stock.cig
  for m in 5 10 20 40 80; do
    local t=()
    for i in 1 2 3; do /usr/bin/time -f "%e" bash -c "MM2_GEO=1 MM2_GEO_MARGIN=$m $MG -ax $preset -t32 $AR $reads >/dev/null 2>/dev/null" 2>$D/t.t; t+=("$(cat $D/t.t)"); done
    MM2_GEO=1 MM2_GEO_MARGIN=$m $MG -ax $preset -t32 $AR $reads 2>/dev/null | cigs > $D/geo.cig
    local cid=$(join -t$'\t' $D/stock.cig $D/geo.cig | awk -F'\t' '{n++;if($2==$3)id++}END{printf "%.2f",(n?100*id/n:0)}')
    local f1="(skip)"
    if [ "$m" = 5 ] || [ "$m" = 20 ] || [ "$m" = 80 ]; then
      MM2_GEO=1 MM2_GEO_MARGIN=$m $MG -ax $preset -t32 $AR $reads 2>/dev/null | q "$SAMTOOLS" sort -@8 -o $D/b.bam -
      q "$SAMTOOLS" reheader -c "sed \"s/SN:chr22/SN:$CT/\"" $D/b.bam > $D/b2.bam; mv $D/b2.bam $D/b.bam; q "$SAMTOOLS" index $D/b.bam
      conda_on "$CONDA_ENV_CLAIR3"; rm -rf $D/c3
      $C3 --bam_fn=$D/b.bam --ref_fn=$CR --threads=32 --platform=$tech --model_path=$mdl --output=$D/c3 --ctg_name=$CT --bed_fn=$BD >$D/c3.log 2>&1
      conda_off; conda_on "$CONDA_ENV_RTG"; rm -rf $D/ev
      if rtg vcfeval -b $TR -c $D/c3/merge_output.vcf.gz -t $SF -e $BD -o $D/ev >/dev/null 2>&1; then
        read ts ti < <(tc $D/ev/tp.vcf.gz); read fps fpi < <(tc $D/ev/fp.vcf.gz); read fns fni < <(tc $D/ev/fn.vcf.gz); read tbs tbi < <(tc $D/ev/tp-baseline.vcf.gz)
        f1=$(awk -v ts=$ts -v ti=$ti -v fps=$fps -v fpi=$fpi -v fns=$fns -v fni=$fni -v tbs=$tbs -v tbi=$tbi 'function F(tp,fp,tb,fn, p,r){p=(tp+fp)?tp/(tp+fp):0;r=(tb+fn)?tb/(tb+fn):0;return (p+r)?2*p*r/(p+r):0} BEGIN{printf "%.4f/%.4f", F(ts,fps,tbs,fns), F(ti,fpi,tbi,fni)}')
      fi
      conda_off
    fi
    awk -v t=$tech -v m=$m -v tm="$(med "${t[@]}")" -v c=$cid -v f="$f1" 'BEGIN{printf "   %-4s m=%-3s time %ss  CIGARid %s%%  SNV/INDEL-F1 %s\n",t,m,tm,c,f}' | tee -a $R
  done
}
mrun hifi map-hifi hifi_chr22.fq $M_HIFI
mrun ont  map-ont  ont_r10_chr22.fq $M_ONT
rm -rf $D; echo "### DONE-MARGIN" | tee -a $R
