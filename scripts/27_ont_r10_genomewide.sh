#!/bin/bash
# ONT R10.4.1 genome-wide subset (Supplementary Tables S21/S22): 10% random sample (samtools -s 42.10, primary
# records) of GIAB HG002 PAO83395.pass.cram, decoded with the GRCh38 no-alt analysis set and aligned to the full
# GRCh38 (map-ont index). Certified identity, 3 interleaved timing runs, DP-stage decomposition.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/ont_r10_genomewide; mkdir -p $D
[ -s "$HG002_ONT_CRAM" ] || { echo "run: bash scripts/fetch_data.sh genomewide"; exit 1; }
[ -s grch38_ont.mmi ] || MM2_GEO=0 "$MM2GEO" -x map-ont -d grch38_ont.mmi "$GRCH38_REF" > /dev/null 2>&1
FQ=ont_r10_genomewide_10pct.fq
[ -s $FQ ] || "$SAMTOOLS" view -@16 -b -T GRCh38_no_alt_analysis_set.fna -F 0x900 -s 42.10 "$HG002_ONT_CRAM" | "$SAMTOOLS" fastq -@8 - > $FQ 2>/dev/null
awk 'NR%4==2{n++;b+=length($0)}END{print "reads "n", bases "b}' $FQ | tee $D/summary.txt
MM2_GEO=0 "$MM2GEO" -cx map-ont -t"$THREADS" grch38_ont.mmi $FQ 2>/dev/null > $D/stock.paf
MM2_GEO=1 MM2_GEO_CERT=1 "$MM2GEO" -cx map-ont -t"$THREADS" grch38_ont.mmi $FQ 2>/dev/null > $D/certified.paf
# differing PAF lines (with CIGAR) and, per differing record, the alignment score (AS) of minimap2 and certified
diff $D/stock.paf $D/certified.paf | grep '^[<>]' > $D/discordant_lines.txt
echo "certified identical to minimap2: $([ -s $D/discordant_lines.txt ] && echo "no, $(cut -c3- $D/discordant_lines.txt | cut -f1 | sort -u | wc -l) reads differ" || echo yes)" | tee -a $D/summary.txt
awk '{for(i=14;i<=NF;i++) if($i~/^AS:i:/) print $1, $2, $7":"$9"-"$10, $i}' $D/discordant_lines.txt | tee -a $D/summary.txt
for rep in 1 2 3; do for c in stock heuristic certified; do
  case $c in stock) E="MM2_GEO=0";; heuristic) E="MM2_GEO=1";; certified) E="MM2_GEO=1 MM2_GEO_CERT=1";; esac
  /usr/bin/time -f "%e %M" -o $D/t env $E "$MM2GEO" -ax map-ont -t"$THREADS" grch38_ont.mmi $FQ > /dev/null 2>/dev/null
  echo -e "$c\trep$rep\t$(cat $D/t)" | tee -a $D/timing.tsv; done; done
for c in stock heuristic certified; do
  case $c in stock) E="MM2_GEO=0";; heuristic) E="MM2_GEO=1";; certified) E="MM2_GEO=1 MM2_GEO_CERT=1";; esac
  /usr/bin/time -f "%e %U %S" -o $D/t env $E MM2_GEO_TIMEDP=1 "$MM2GEO_INSTR" -ax map-ont -t"$THREADS" grch38_ont.mmi $FQ > /dev/null 2> $D/err
  read e u s < $D/t; echo -e "$c\twall=$e\tcpu=$(echo "$u + $s" | bc)\tdp_cpu=$(grep -o 'dp_cpu_s=[0-9.]*' $D/err | cut -d= -f2)" | tee -a $D/dpstage.tsv; done
