#!/bin/bash
# Large-SV band-stress: does geo's speedup shrink toward 1x (band -> cap) while staying
# alignment-concordant with stock? L=300-1500bp dense-ish. clair3 can't score big indels,
# so accuracy proxy = CIGAR identity of geo vs stock on shared primary alignments.
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"
BGZIP="$BGZIP"
TABIX="$TABIX"
BCFTOOLS="$BCFTOOLS"
source "$CONDA_PROFILE"
export LC_ALL=C; q(){ "$@" 2>/dev/null; }
D=advsv_chr22; mkdir -p $D
R=hifi_bench/results_ADVSV.txt; : > $R
med(){ printf '%s\n' "$@"|sort -n|awk '{a[NR]=$1}END{print a[int((NR+1)/2)]}'; }

echo ">> [1] plant LARGE indels (300-1500bp)" | tee -a $R
ADV_OUT=$D/truth_sv.vcf ADV_STEP=3000 ADV_LMIN=300 ADV_LMAX=1500 ADV_SEED=77 python3 make_adversarial.py | tee -a $R
$BGZIP -f $D/truth_sv.vcf; $TABIX -f -p vcf $D/truth_sv.vcf.gz

echo ">> [2] build large-SV sample (bcftools consensus)" | tee -a $R
$BCFTOOLS consensus -f chr22_named.fa $D/truth_sv.vcf.gz > $D/sample_sv.fa 2>>$R
echo "   sample bp: $(grep -v '^>' $D/sample_sv.fa | tr -d '\n' | wc -c)" | tee -a $R

echo ">> [3] simulate HiFi 15x" | tee -a $R
badread simulate --seed 77 --reference $D/sample_sv.fa --quantity 15x --length 12000,6000 \
   --identity 99.5,100,1 --junk_reads 0 --random_reads 0 --chimeras 0 > $D/sv.fq 2>/dev/null
echo "   reads: $(($(wc -l < $D/sv.fq)/4))" | tee -a $R

echo ">> [4] alignment time (index build included), median of 3 (32t)" | tee -a $R
for cfg in mm2 geo; do
  gf=0; [ "$cfg" = geo ] && gf=1; t=()
  for i in 1 2 3; do
    /usr/bin/time -f "%e" bash -c "MM2_GEO=$gf MM2_GEO_MARGIN=20 $MG -ax map-hifi -t32 chr22_named.fa $D/sv.fq >/dev/null 2>/dev/null" 2>$D/t.t
    t+=("$(cat $D/t.t)")
  done
  printf "   %-4s time median: %ss  (runs: %s)\n" "$cfg" "$(med "${t[@]}")" "${t[*]}" | tee -a $R
done

echo ">> [5] CIGAR concordance geo vs stock (accuracy proxy; clair3 can't score big indels)" | tee -a $R
# primary alignments only (-F 0x900), name-sorted, compare CIGAR per shared read
MM2_GEO=0 $MG -ax map-hifi -t32 chr22_named.fa $D/sv.fq 2>/dev/null | q "$SAMTOOLS" view -F 0x900 - | awk '{print $1"\t"$6}' | sort > $D/mm2.cig
MM2_GEO=1 MM2_GEO_MARGIN=20 $MG -ax map-hifi -t32 chr22_named.fa $D/sv.fq 2>/dev/null | q "$SAMTOOLS" view -F 0x900 - | awk '{print $1"\t"$6}' | sort > $D/geo.cig
join -t$'\t' $D/mm2.cig $D/geo.cig | awk -F'\t' '{n++; if($2==$3)id++} END{printf "   shared primaries: %d  CIGAR-identical: %d (%.2f%%)\n", n, id, (n?100*id/n:0)}' | tee -a $R
echo "### DONE-ADVSV" | tee -a $R
