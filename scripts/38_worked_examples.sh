#!/bin/bash
# Worked examples of the heuristic mode's limitation (Supplementary Section "Localisation of the CIGAR-discordant
# reads", paragraph "Two examples"): six chr22 HiFi reads aligned with minimap2 (MM2_GEO=0), the heuristic mode
# (MM2_GEO=1 MM2_GEO_FAST=1) and the certified mode (MM2_GEO=1), one thread, PAF with CIGAR; plus the per-gap log
# (MM2_GEO_GAPLOG, instrumented build) of read m64011_190901_095311/1378531/ccs in the three configurations.
# Output: results/examples/{reads.txt,examples_{stock,heur,cert}.paf,gaplog_1378531_{stock,heur,cert}.tsv}
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/examples; mkdir -p $D
[ -s chr22.hifi.mmi ] || MM2_GEO=0 "$MM2GEO" -x map-hifi -d chr22.hifi.mmi chr22_named.fa > /dev/null 2>&1
cat > $D/reads.txt <<'EOF'
m64011_190830_220126/18024473/ccs
m64015_190922_010918/131335307/ccs
m64011_190830_220126/34801958/ccs
m64011_190830_220126/76350446/ccs
m64011_190901_095311/1378531/ccs
m64015_190922_010918/68879252/ccs
EOF
pick(){ awk 'NR==FNR{w["@"$1]=1;next} FNR%4==1{k=($1 in w)} k' "$1" hifi_chr22.fq; }   # names file -> FASTQ records
pick $D/reads.txt > $D/ex.fq
echo m64011_190901_095311/1378531/ccs > $D/r1.txt; pick $D/r1.txt > $D/r1.fq
for m in stock heur cert; do
  case $m in stock) E="MM2_GEO=0";; heur) E="MM2_GEO=1 MM2_GEO_FAST=1";; cert) E="MM2_GEO=1";; esac
  env $E "$MM2GEO" -cx map-hifi -t1 chr22.hifi.mmi $D/ex.fq > $D/examples_$m.paf 2>/dev/null
  env $E MM2_GEO_GAPLOG=$D/gaplog_1378531_$m.tsv "$MM2GEO_INSTR" -cx map-hifi -t1 chr22.hifi.mmi $D/r1.fq > /dev/null 2>&1
done
rm -f $D/ex.fq $D/r1.fq $D/r1.txt
for m in stock heur cert; do   # primary alignments: read, query span, strand:target start, AS, CIGAR (truncated)
  awk -v m=$m '/tp:A:P/{for(i=13;i<=NF;i++){if($i~/^AS:/)a=$i; if($i~/^cg:/)c=$i}; print m, $1, $3"-"$4, $5":"$8, a, substr(c,1,90)}' $D/examples_$m.paf
done
