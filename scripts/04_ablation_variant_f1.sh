#!/bin/bash
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"
C3="$CLAIR3"
M_HIFI="$MODEL_HIFI"
M_ONT="$MODEL_ONT"
AR=chr22_named.fa; CR=chr22_cm.fa; CT=CM000684.2; TR=truth_chr22.vcf.gz; BD=conf_chr22.bed; SF=ref_chr22_sdf
source "$CONDA_PROFILE"
export LC_ALL=C; q(){ "$@" 2>/dev/null; }
D=ablf1; mkdir -p $D; R=hifi_bench/results_ABLATION_F1.txt; : > $R
tc(){ zcat "$1" 2>/dev/null | awk '!/^#/{if(length($4)==1&&length($5)==1)s++;else i++}END{printf "%d %d",s+0,i+0}'; }
f1(){ local tag=$1 tech=$2 preset=$3 reads=$4 mdl=$5 env=$6
  eval "MM2_GEO=1 $env $MG -ax $preset -t32 $AR $reads 2>/dev/null" | q "$SAMTOOLS" sort -@8 -o $D/b.bam -
  q "$SAMTOOLS" reheader -c "sed \"s/SN:chr22/SN:$CT/\"" $D/b.bam > $D/b2.bam; mv $D/b2.bam $D/b.bam; q "$SAMTOOLS" index $D/b.bam
  conda_on "$CONDA_ENV_CLAIR3"; rm -rf $D/c3
  $C3 --bam_fn=$D/b.bam --ref_fn=$CR --threads=32 --platform=$tech --model_path=$mdl --output=$D/c3 --ctg_name=$CT --bed_fn=$BD >$D/c3.log 2>&1
  conda_off; conda_on "$CONDA_ENV_RTG"; rm -rf $D/ev
  rtg vcfeval -b $TR -c $D/c3/merge_output.vcf.gz -t $SF -e $BD -o $D/ev >/dev/null 2>&1 || { echo "$tag: rtg vcfeval FAILED (see $D/c3.log)" | tee -a $R; conda_off; return; }
  read ts ti < <(tc $D/ev/tp.vcf.gz); read fps fpi < <(tc $D/ev/fp.vcf.gz); read fns fni < <(tc $D/ev/fn.vcf.gz); read tbs tbi < <(tc $D/ev/tp-baseline.vcf.gz)
  conda_off
  awk -v t="$tag" -v ts=$ts -v ti=$ti -v fps=$fps -v fpi=$fpi -v fns=$fns -v fni=$fni -v tbs=$tbs -v tbi=$tbi 'function F(tp,fp,tb,fn,p,r){p=(tp+fp)?tp/(tp+fp):0;r=(tb+fn)?tb/(tb+fn):0;return (p+r)?2*p*r/(p+r):0}BEGIN{printf "%-40s SNV/INDEL %.4f/%.4f\n",t,F(ts,fps,tbs,fns),F(ti,fpi,tbi,fni)}' | tee -a $R
}
for TE in "hifi map-hifi hifi_chr22.fq $M_HIFI" "ont map-ont ont_r10_chr22.fq $M_ONT"; do
  set -- $TE; tech=$1; pre=$2; rd=$3; mdl=$4
  echo "==== $tech ====" | tee -a $R
  f1 "$tech fixed-narrow"   $tech $pre $rd $mdl "MM2_GEO_FIXED=1 MM2_GEO_NODOUBLE=1"
  f1 "$tech fixed-doubling" $tech $pre $rd $mdl "MM2_GEO_FIXED=1"
  f1 "$tech geometry-only"  $tech $pre $rd $mdl "MM2_GEO_NODOUBLE=1"
  f1 "$tech mm2-geo"        $tech $pre $rd $mdl ""
done
rm -rf $D; echo "### DONE-ABLF1" | tee -a $R
