#!/bin/bash
# Three-caller cross-check (Supplementary Table S5): HiFi chr22 aligned by minimap2 and by the heuristic mode
# (m=20), small variants called with Clair3 (local re-alignment), bcftools mpileup+call (pileup-based, trusts
# the CIGAR) and DeepVariant 1.6.1 (docker image google/deepvariant:1.6.1), each scored with rtg vcfeval against
# GIAB v4.2.1 in the high-confidence regions. Needs docker for DeepVariant (set DEEPVARIANT=skip to omit it).
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"
DV_IMAGE="${DV_IMAGE:-google/deepvariant:1.6.1}"; DEEPVARIANT="${DEEPVARIANT:-docker}"
S=tc_scratch; mkdir -p $S; R=hifi_bench/results_THREE_CALLERS.txt; : > $R
export LC_ALL=C
tc(){ zcat "$1" 2>/dev/null | awk '!/^#/{if(length($4)==1&&length($5)==1)s++;else i++}END{printf "%d %d",s+0,i+0}'; }
score(){ # label vcf
  conda_on "$CONDA_ENV_RTG"; rm -rf $S/ev
  "$RTG" vcfeval -b truth_chr22.vcf.gz -c "$2" -t ref_chr22_sdf -e conf_chr22.bed -o $S/ev >/dev/null 2>&1 || { echo "  $1 FAIL" | tee -a $R; conda_off; return; }
  read ts ti < <(tc $S/ev/tp.vcf.gz); read fps fpi < <(tc $S/ev/fp.vcf.gz); read fns fni < <(tc $S/ev/fn.vcf.gz); read tbs tbi < <(tc $S/ev/tp-baseline.vcf.gz)
  conda_off
  awk -v ts=$ts -v ti=$ti -v fps=$fps -v fpi=$fpi -v fns=$fns -v fni=$fni -v tbs=$tbs -v tbi=$tbi -v l="$1" \
   'function f(tp,fp,tb,fn, p,r){p=(tp+fp)?tp/(tp+fp):0;r=(tb+fn)?tb/(tb+fn):0;return sprintf("%.4f/%.4f/%.4f",p,r,(p+r)?2*p*r/(p+r):0)}
    BEGIN{printf "  %-24s SNV P/R/F1 %-22s INDEL P/R/F1 %-22s\n",l,f(ts,fps,tbs,fns),f(ti,fpi,tbi,fni)}' | tee -a $R; }
echo "### HiFi chr22: minimap2 vs heuristic mode (m=20), three callers" | tee -a $R
for cfg in mm2 geo; do
  E="MM2_GEO=0"; [ $cfg = geo ] && E="MM2_GEO=1 MM2_GEO_MARGIN=20"
  env $E "$MM2GEO" -ax map-hifi -t"$THREADS" chr22_named.fa hifi_chr22.fq 2>/dev/null | "$SAMTOOLS" sort -@8 -o $S/tmp.bam - 2>/dev/null
  "$SAMTOOLS" reheader -c 'sed "s/SN:chr22/SN:CM000684.2/"' $S/tmp.bam > $S/$cfg.bam 2>/dev/null; rm -f $S/tmp.bam
  "$SAMTOOLS" index -@8 $S/$cfg.bam
  # Clair3
  conda_on "$CONDA_ENV_CLAIR3"; rm -rf $S/c3_$cfg
  "$CLAIR3" --bam_fn=$S/$cfg.bam --ref_fn=chr22_cm.fa --threads="$THREADS" --platform=hifi --model_path="$MODEL_HIFI" \
     --output=$S/c3_$cfg --ctg_name=CM000684.2 --bed_fn=conf_chr22.bed > $S/c3_$cfg.log 2>&1
  conda_off; score "clair3 $cfg" $S/c3_$cfg/merge_output.vcf.gz
  # bcftools mpileup + call
  "$BCFTOOLS" mpileup -f chr22_cm.fa -R conf_chr22.bed -a AD,DP -Ou $S/$cfg.bam 2>/dev/null | "$BCFTOOLS" call -mv -Oz -o $S/bcf_$cfg.vcf.gz 2>/dev/null
  "$TABIX" -f -p vcf $S/bcf_$cfg.vcf.gz; score "bcftools $cfg" $S/bcf_$cfg.vcf.gz
  # DeepVariant
  if [ "$DEEPVARIANT" = docker ]; then
    docker run --rm -u "$(id -u):$(id -g)" -v "$WORKDIR":/data "$DV_IMAGE" /opt/deepvariant/bin/run_deepvariant \
      --model_type=PACBIO --ref=/data/chr22_cm.fa --reads=/data/$S/$cfg.bam --regions=/data/conf_chr22.bed \
      --output_vcf=/data/$S/dv_$cfg.vcf.gz --num_shards="$THREADS" --intermediate_results_dir=/data/$S/dvtmp_$cfg > $S/dv_$cfg.log 2>&1
    score "deepvariant $cfg" $S/dv_$cfg.vcf.gz
  fi
done
rm -rf $S; echo "### DONE" | tee -a $R
