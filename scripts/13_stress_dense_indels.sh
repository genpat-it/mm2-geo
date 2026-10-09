#!/bin/bash
# Adversarial indel-dense stress test on chr22: does geo's speedup collapse (graceful) while F1 holds?
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"
BGZIP="$BGZIP"
TABIX="$TABIX"
BCFTOOLS="$BCFTOOLS"
C3="$CLAIR3"
M_HIFI="$MODEL_HIFI"
source "$CONDA_PROFILE"
export LC_ALL=C; q(){ "$@" 2>/dev/null; }
D=adv_chr22; mkdir -p $D
R=hifi_bench/results_ADVERSARIAL.txt; : > $R
med(){ printf '%s\n' "$@"|sort -n|awk '{a[NR]=$1}END{print a[int((NR+1)/2)]}'; }

echo ">> [1] plant dense indels" | tee -a $R
ADV_SEED=77 python3 "$SCRIPTS/make_adversarial.py" | tee -a $R
$BGZIP -f $D/truth_adv.vcf; $TABIX -f -p vcf $D/truth_adv.vcf.gz

echo ">> [2] build adversarial sample genome (bcftools consensus)" | tee -a $R
$BCFTOOLS consensus -f chr22_named.fa $D/truth_adv.vcf.gz > $D/sample_adv.fa 2>>$R
echo "   sample bp: $(grep -v '^>' $D/sample_adv.fa | tr -d '\n' | wc -c)" | tee -a $R

echo ">> [3] simulate HiFi 20x (badread)" | tee -a $R
badread simulate --seed 77 --reference $D/sample_adv.fa --quantity 20x --length 12000,6000 \
   --identity 99.5,100,1 --junk_reads 0 --random_reads 0 --chimeras 0 > $D/adv.fq 2>/dev/null
echo "   reads: $(($(wc -l < $D/adv.fq)/4))" | tee -a $R

echo ">> [4] index ref + build SDF" | tee -a $R
q "$SAMTOOLS" faidx chr22_named.fa
conda_on "$CONDA_ENV_RTG"; rm -rf $D/sdf; rtg format -o $D/sdf chr22_named.fa >/dev/null 2>&1; conda_off

echo ">> [5] alignment time (pre-built map-hifi index, so index construction is not timed), median of 3 (32t)" | tee -a $R
[ -s chr22.hifi.mmi ] || $MG -x map-hifi -d chr22.hifi.mmi chr22_named.fa >/dev/null 2>&1
for cfg in mm2 geo cert; do
  case $cfg in mm2) E="MM2_GEO=0";; geo) E="MM2_GEO=1 MM2_GEO_MARGIN=20";; cert) E="MM2_GEO=1 MM2_GEO_CERT=1";; esac
  t=()
  for i in 1 2 3; do
    /usr/bin/time -f "%e" bash -c "$E $MG -ax map-hifi -t32 chr22.hifi.mmi $D/adv.fq >/dev/null 2>/dev/null" 2>$D/t.t
    t+=("$(cat $D/t.t)")
  done
  printf "   %-4s time median: %ss  (runs: %s)\n" "$cfg" "$(med "${t[@]}")" "${t[*]}" | tee -a $R
done

echo ">> [6] F1 (clair3 hifi + rtg vcfeval vs planted indels)" | tee -a $R
tc(){ zcat "$1" 2>/dev/null | awk '!/^#/{if(length($4)==1&&length($5)==1)s++;else i++}END{printf "%d %d",s+0,i+0}'; }
for cfg in mm2 geo cert; do
  case $cfg in mm2) E="MM2_GEO=0";; geo) E="MM2_GEO=1 MM2_GEO_MARGIN=20";; cert) E="MM2_GEO=1 MM2_GEO_CERT=1";; esac
  env $E $MG -ax map-hifi -t32 chr22_named.fa $D/adv.fq 2>/dev/null | q "$SAMTOOLS" sort -@8 -o $D/b.bam -
  q "$SAMTOOLS" index $D/b.bam
  conda_on "$CONDA_ENV_CLAIR3"; rm -rf $D/c3
  $C3 --bam_fn=$D/b.bam --ref_fn=chr22_named.fa --threads=32 --platform=hifi --model_path=$M_HIFI --output=$D/c3 --ctg_name=chr22 --include_all_ctgs --enable_long_indel >$D/c3.log 2>&1
  conda_off; conda_on "$CONDA_ENV_RTG"; rm -rf $D/ev
  if ! rtg vcfeval -b $D/truth_adv.vcf.gz -c $D/c3/merge_output.vcf.gz -t $D/sdf --squash-ploidy -o $D/ev >/dev/null 2>&1; then
    echo "   $cfg vcfeval FAIL" | tee -a $R; conda_off; continue; fi
  read ts ti < <(tc $D/ev/tp.vcf.gz); read fps fpi < <(tc $D/ev/fp.vcf.gz); read fns fni < <(tc $D/ev/fn.vcf.gz); read tbs tbi < <(tc $D/ev/tp-baseline.vcf.gz)
  conda_off
  awk -v ts=$ts -v ti=$ti -v fps=$fps -v fpi=$fpi -v fns=$fns -v fni=$fni -v tbs=$tbs -v tbi=$tbi -v c="$cfg" \
   'function f(tp,fp,tb,fn, p,r){p=(tp+fp)?tp/(tp+fp):0;r=(tb+fn)?tb/(tb+fn):0;return sprintf("%.4f/%.4f/%.4f",p,r,(p+r)?2*p*r/(p+r):0)}
    BEGIN{printf "   %-4s SNV %-22s INDEL %-22s\n",c,f(ts,fps,tbs,fns),f(ti,fpi,tbi,fni)}' | tee -a $R
done
echo "### DONE-ADV" | tee -a $R
