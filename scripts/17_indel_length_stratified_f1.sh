#!/bin/bash
# Supplementary Table S10: small-variant F1 stratified by indel length.
# Stock minimap2 (MM2_GEO=0) vs mm2-geo (MM2_GEO=1) on the SAME binary, HiFi + ONT chr22.
# Runs Clair3 + rtg vcfeval retaining tp-baseline/fp/fn, then 17_stratify_indel.py bins by
# indel length L=|len(REF)-len(ALT)| into SNV / 1 / 2-5 / 6-15 / 16-50 bp and recomputes F1.
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"                       # instrumented build == ./mm2geo (see tool repo)
C3="$CLAIR3"
M_HIFI="$MODEL_HIFI"
M_ONT="$MODEL_ONT"
AR=chr22_named.fa; CR=chr22_cm.fa; CT=CM000684.2
TR=truth_chr22.vcf.gz; BD=conf_chr22.bed; SF=ref_chr22_sdf
export TMPDIR="$WORKDIR"/strat_tmp
mkdir -p "$TMPDIR"
source "$CONDA_PROFILE"
export LC_ALL=C; q(){ "$@" 2>/dev/null; }
OUT=strat; mkdir -p $OUT
R=hifi_bench/results_STRAT_INDEL.txt; : > $R

run(){ local tech=$1 preset=$2 reads=$3 mdl=$4 aln=$5 geo=$6   # aln: mm2|geo ; geo: 0|1
  local tag="${tech}_${aln}"; local W=$OUT/$tag
  echo ">>> $tag (MM2_GEO=$geo)" | tee -a $R
  rm -rf $W; mkdir -p $W
  MM2_GEO=$geo $MG -ax $preset -t32 $AR $reads 2>/dev/null | q "$SAMTOOLS" sort -@8 -T $TMPDIR/srt_$tag -o $W/b.bam -
  q "$SAMTOOLS" reheader -c "sed \"s/SN:chr22/SN:$CT/\"" $W/b.bam > $W/b2.bam; mv $W/b2.bam $W/b.bam
  q "$SAMTOOLS" index $W/b.bam
  conda activate "$CONDA_ENV_CLAIR3"; rm -rf $W/c3
  $C3 --bam_fn=$W/b.bam --ref_fn=$CR --threads=32 --platform=$tech --model_path=$mdl \
      --output=$W/c3 --ctg_name=$CT --bed_fn=$BD >$W/c3.log 2>&1
  conda deactivate; conda activate "$CONDA_ENV_RTG"; rm -rf $W/ev
  rtg vcfeval -b $TR -c $W/c3/merge_output.vcf.gz -t $SF -e $BD -o $W/ev >/dev/null 2>&1 || { echo "    $tag: rtg vcfeval FAILED (see $W/c3.log)" | tee -a $R; conda deactivate; return; }
  conda deactivate
  echo "    eval retained: $W/ev (tp-baseline/fp/fn)" | tee -a $R
}

for TE in "hifi map-hifi hifi_chr22.fq $M_HIFI" "ont map-ont ont_r10_chr22.fq $M_ONT"; do
  set -- $TE; tech=$1; pre=$2; rd=$3; mdl=$4
  run $tech $pre $rd $mdl mm2 0
  run $tech $pre $rd $mdl geo 1
done

echo "### DONE-STRAT" | tee -a $R
echo | tee -a $R
python3 "$SCRIPTS/17_stratify_indel.py" strat | tee -a $R
