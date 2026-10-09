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
ADV_OUT=$D/truth_sv.vcf ADV_STEP=3000 ADV_LMIN=300 ADV_LMAX=1500 ADV_SEED=77 python3 "$SCRIPTS/make_adversarial.py" | tee -a $R
$BGZIP -f $D/truth_sv.vcf; $TABIX -f -p vcf $D/truth_sv.vcf.gz

echo ">> [2] build large-SV sample (bcftools consensus)" | tee -a $R
$BCFTOOLS consensus -f chr22_named.fa $D/truth_sv.vcf.gz > $D/sample_sv.fa 2>>$R
echo "   sample bp: $(grep -v '^>' $D/sample_sv.fa | tr -d '\n' | wc -c)" | tee -a $R

echo ">> [3] simulate HiFi 15x" | tee -a $R
badread simulate --seed 77 --reference $D/sample_sv.fa --quantity 15x --length 12000,6000 \
   --identity 99.5,100,1 --junk_reads 0 --random_reads 0 --chimeras 0 > $D/sv.fq 2>/dev/null
echo "   reads: $(($(wc -l < $D/sv.fq)/4))" | tee -a $R

echo ">> [4] alignment time (pre-built map-hifi index, so index construction is not timed), median of 3 (32t)" | tee -a $R
[ -s chr22.hifi.mmi ] || $MG -x map-hifi -d chr22.hifi.mmi chr22_named.fa >/dev/null 2>&1
env_of(){ case $1 in mm2) echo "MM2_GEO=0";; geo) echo "MM2_GEO=1 MM2_GEO_MARGIN=20";; cert) echo "MM2_GEO=1 MM2_GEO_CERT=1";; esac; }
for cfg in mm2 geo cert; do
  t=()
  for i in 1 2 3; do
    /usr/bin/time -f "%e" bash -c "$(env_of $cfg) $MG -ax map-hifi -t32 chr22.hifi.mmi $D/sv.fq >/dev/null 2>/dev/null" 2>$D/t.t
    t+=("$(cat $D/t.t)")
  done
  printf "   %-4s time median: %ss  (runs: %s)\n" "$cfg" "$(med "${t[@]}")" "${t[*]}" | tee -a $R
done

echo ">> [5] placement and CIGAR concordance vs stock, heuristic (geo) and certified (cert) (accuracy proxy; clair3 can't score big indels)" | tee -a $R
# primary alignments only (-F 0x900): read, mapped flag, contig, strand, position, CIGAR; name-sorted
ex(){ q "$SAMTOOLS" view -F 0x900 - | awk 'BEGIN{OFS="\t"}{print $1, (and($2,4)?"u":"m"), $3, (and($2,16)?"-":"+"), $4, $6}' | LC_ALL=C sort -t$'\t' -k1,1; }
for cfg in mm2 geo cert; do env $(env_of $cfg) $MG -ax map-hifi -t32 chr22.hifi.mmi $D/sv.fq 2>/dev/null | ex > $D/$cfg.prim; done
for cfg in geo cert; do
  LC_ALL=C join -t$'\t' $D/mm2.prim $D/$cfg.prim | awk -F'\t' -v c=$cfg '{n++; if($2=="m")ma++; if($7=="m")mb++; if($2=="m" && $7=="m"){both++; if($3==$8 && $4==$9 && $5==$10)pl++; if($6==$11)id++}}
    END{printf "   %-4s reads %d  mapped mm2/%s %d/%d  of reads mapped by both (%d): same placement (contig+strand+POS) %.2f%%, CIGAR-identical %.2f%%\n", c, n, c, ma, mb, both, 100*pl/both, 100*id/both}' | tee -a $R
done
echo "### DONE-ADVSV" | tee -a $R
