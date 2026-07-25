#!/bin/bash
# Table 7 / Supplementary Table S8: structural-variant F1, stock minimap2 (MM2_GEO=0) vs mm2-geo
# (MM2_GEO=1) on the SAME binary -> Sniffles2 -> truvari, vs GIAB HG002 SV Tier1 v0.6 (GRCh37 chr22).
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
S="$WORKDIR/sv"; mkdir -p "$S"
R="$WORKDIR/hifi_bench/results_SV_geo.txt"; mkdir -p "$WORKDIR/hifi_bench"; : > "$R"

# 1) reference: GRCh37 chr22 (Ensembl). GIAB SV Tier1 truth uses GRCh37.
if [ ! -s "$S/grch37_chr22.fa" ]; then
  echo ">> downloading GRCh37 chr22 (Ensembl)" | tee -a "$R"
  curl -sSL "http://ftp.ensembl.org/pub/grch37/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz" \
    | gunzip -c > "$S/grch37_chr22.fa" || { echo "GRCh37 download FAILED" | tee -a "$R"; exit 1; }
fi
"$SAMTOOLS" faidx "$S/grch37_chr22.fa"
REF="$S/grch37_chr22.fa"
CTG=$(head -1 "$REF" | sed 's/>//;s/ .*//')

# 2) restrict GIAB HG002 SV Tier1 v0.6 truth + high-confidence BED to this contig.
#    Place HG002_SVs_Tier1_v0.6.vcf.gz and .bed in $S first (see ACCESSIONS.txt).
TR="$S/truth_${CTG}.vcf.gz"; BED="$S/tier1_${CTG}.bed"
"$BCFTOOLS" view "$S/HG002_SVs_Tier1_v0.6.vcf.gz" -r "$CTG" -Oz -o "$TR" 2>/dev/null \
  || "$BCFTOOLS" view "$S/HG002_SVs_Tier1_v0.6.vcf.gz" -Oz -o "$TR"
"$BCFTOOLS" index -tf "$TR"
awk -v c="$CTG" '$1==c' "$S/HG002_SVs_Tier1_v0.6.bed" > "$BED"

# 3) per technology: map with stock (gf=0) and geo (gf=1), call SVs, benchmark.
run(){ local tech=$1 preset=$2 reads=$3
  echo ">> $tech (SV F1 vs GIAB Tier1, GRCh37 chr22)" | tee -a "$R"
  for cfg in mm2 geo; do
    local gf=0; [ "$cfg" = geo ] && gf=1
    MM2_GEO=$gf MM2_GEO_MARGIN=20 "$MM2GEO" -ax "$preset" -t"$THREADS" "$REF" "$WORKDIR/$reads" 2>/dev/null \
      | "$SAMTOOLS" sort -@8 -o "$S/${tech}_${cfg}.bam" - 2>/dev/null
    "$SAMTOOLS" index "$S/${tech}_${cfg}.bam"
    "$SNIFFLES" --input "$S/${tech}_${cfg}.bam" --vcf "$S/${tech}_${cfg}_sv.vcf.gz" --reference "$REF" --threads 16 \
      >"$S/snf_${tech}_${cfg}.log" 2>&1
    echo "   $cfg raw SV calls: $(zcat "$S/${tech}_${cfg}_sv.vcf.gz" 2>/dev/null | grep -vc '^#')" | tee -a "$R"
    rm -rf "$S/${tech}_${cfg}_truvari"
    "$TRUVARI" bench -b "$TR" -c "$S/${tech}_${cfg}_sv.vcf.gz" --includebed "$BED" -o "$S/${tech}_${cfg}_truvari" \
      >"$S/truv_${tech}_${cfg}.log" 2>&1 || true
    local j="$S/${tech}_${cfg}_truvari/summary.json"
    if [ -f "$j" ]; then
      python3 -c "import json;d=json.load(open('$j'));print('   %-4s %-3s  P=%.3f R=%.3f F1=%.3f  TP=%d FP=%d FN=%d'%('$tech','$cfg',d.get('precision',0),d.get('recall',0),d.get('f1',0),d.get('TP-comp',d.get('TP',0)),d.get('FP',0),d.get('FN',0)))" | tee -a "$R"
    else echo "   $tech $cfg truvari FAIL (see truv_${tech}_${cfg}.log)" | tee -a "$R"; fi
  done
}
run ONT  map-ont  ont_r10_chr22.fq
run HiFi map-hifi hifi_chr22.fq
echo "### DONE-SVGEO" | tee -a "$R"
