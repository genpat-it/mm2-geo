#!/bin/bash
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"; SAM="$SAMTOOLS"
C3="$CLAIR3"
MODEL="$MODEL_BACT"
REF=ecoli_vc/ref_mut.fa; TR=ecoli_vc/truth.vcf.gz; RD=ecoli_real2/SRR9900640.fastq
source "$CONDA_PROFILE"
export LC_ALL=C TMPDIR="$WORKDIR"/ecoli_vc/tmp
mkdir -p $TMPDIR
R=hifi_bench/results_ECOLI_F1.txt; : > $R
q(){ "$@" 2>/dev/null; }
$SAM faidx $REF
# rtg SDF
conda activate "$CONDA_ENV_RTG"; rm -rf ecoli_vc/sdf; rtg format -o ecoli_vc/sdf $REF >/dev/null 2>&1; conda deactivate
tc(){ zcat "$1" 2>/dev/null | awk '!/^#/{if(length($4)==1&&length($5)==1)s++;else i++}END{printf "%d %d",s+0,i+0}'; }
f1(){ local tag=$1 G=$2
  MM2_GEO=$G MM2_GEO_MARGIN=20 $MG -ax map-ont -t32 $REF $RD 2>/dev/null | q $SAM sort -@8 -o ecoli_vc/$tag.bam -
  q $SAM index ecoli_vc/$tag.bam
  conda activate "$CONDA_ENV_CLAIR3"; rm -rf ecoli_vc/c3_$tag
  $C3 --bam_fn=ecoli_vc/$tag.bam --ref_fn=$REF --threads=32 --platform=ont \
    --model_path=$MODEL --output=ecoli_vc/c3_$tag --include_all_ctgs --haploid_precise \
    --no_phasing_for_fa >ecoli_vc/c3_$tag.log 2>&1
  conda deactivate; conda activate "$CONDA_ENV_RTG"; rm -rf ecoli_vc/ev_$tag
  rtg vcfeval -b $TR -c ecoli_vc/c3_$tag/merge_output.vcf.gz -t ecoli_vc/sdf -o ecoli_vc/ev_$tag --squash-ploidy >/dev/null 2>&1 || { echo "$tag: rtg vcfeval FAILED (see ecoli_vc/c3_$tag.log)" | tee -a $R; conda deactivate; return; }
  read ts ti < <(tc ecoli_vc/ev_$tag/tp.vcf.gz); read fps fpi < <(tc ecoli_vc/ev_$tag/fp.vcf.gz)
  read fns fni < <(tc ecoli_vc/ev_$tag/fn.vcf.gz); read tbs tbi < <(tc ecoli_vc/ev_$tag/tp-baseline.vcf.gz)
  conda deactivate
  awk -v t=$tag -v ts=$ts -v ti=$ti -v fps=$fps -v fpi=$fpi -v fns=$fns -v fni=$fni -v tbs=$tbs -v tbi=$tbi 'function F(tp,fp,tb,fn,p,r){p=(tp+fp)?tp/(tp+fp):0;r=(tb+fn)?tb/(tb+fn):0;return (p+r)?2*p*r/(p+r):0}BEGIN{printf "%-8s SNV/INDEL %.4f/%.4f  (tpSNP=%d fpSNP=%d fnSNP=%d | tpI=%d fpI=%d fnI=%d)\n",t,F(ts,fps,tbs,fns),F(ti,fpi,tbi,fni),ts,fps,fns,ti,fpi,fni}' | tee -a $R
}
f1 mm2 0
f1 geo 1
echo "### DONE-ECOLIF1" | tee -a $R
