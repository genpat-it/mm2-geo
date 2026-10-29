#!/bin/bash
# ONT R9.4.1 genome-wide subset (Supplementary Tables S19, S21, S22): the first 20 GB of the GIAB HG002 ultra-long
# Guppy 3.2.4 FASTQ (1,228,750 reads, 9.90 Gb; streamed by fetch_data.sh genomewide), aligned to the full GRCh38
# (map-ont index). Certified identity, 3 interleaved timing runs, DP-stage decomposition.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/ont_r9_genomewide; mkdir -p $D
: > $D/timing.tsv; : > $D/dpstage.tsv   # fresh tables on every run (summary.txt is rewritten below)
[ -s "$ONT_R9_READS" ] || { echo "run: bash scripts/fetch_data.sh genomewide"; exit 1; }
[ -s grch38_ont.mmi ] || MM2_GEO=0 "$MM2GEO" -x map-ont -d grch38_ont.mmi "$GRCH38_REF" > /dev/null 2>&1
FQ="$ONT_R9_READS"
awk 'NR%4==2{n++;b+=length($0)}END{print "reads "n", bases "b}' $FQ | tee $D/summary.txt
MM2_GEO=0 "$MM2GEO" -cx map-ont -t"$THREADS" grch38_ont.mmi $FQ 2>/dev/null > $D/stock.paf
# mapped reads (reads with >= 1 primary record, tp:A:P) and primary records (supplementary alignments included,
# which minimap2 also tags tp:A:P): the definition of Supplementary Table S19 (script 29)
awk -v n="$(awk 'NR==1{print $2+0}' $D/summary.txt)" '$0~/tp:A:P/{p++; r[$1]=1} END{m=length(r); printf "mapped reads %d (%.2f%%), primary records %d\n", m, 100*m/n, p}' $D/stock.paf | tee -a $D/summary.txt
# mapped fraction by read length (the R9.4.1 reads are unfiltered; Supplementary Table S19, note a)
awk '$0~/tp:A:P/{print $1}' $D/stock.paf | sort -u > $D/mapped_names.txt
awk 'NR==FNR{m[$1]=1;next} FNR%4==1{n=substr($1,2); sub(/ .*/,"",n)} FNR%4==2{L=length($0); b=(L<1000?"<1k":L<5000?"1-5k":L<20000?"5-20k":">=20k"); t[b]++; if(n in m){c[b]++; MB+=L}; TB+=L}
  END{for(b in t) printf "   reads %s: %d, mapped %.1f%%\n", b, t[b], 100*c[b]/t[b]; printf "   bases in mapped reads %.1f%%\n", 100*MB/TB}' $D/mapped_names.txt $FQ | tee -a $D/summary.txt
rm -f $D/mapped_names.txt
MM2_GEO=1 "$MM2GEO" -cx map-ont -t"$THREADS" grch38_ont.mmi $FQ 2>/dev/null > $D/certified.paf
# differing PAF lines (with CIGAR) and, per differing record, the alignment score (AS) of minimap2 and certified
diff $D/stock.paf $D/certified.paf | grep '^[<>]' > $D/discordant_lines.txt
echo "certified identical to minimap2: $([ -s $D/discordant_lines.txt ] && echo "no, $(cut -c3- $D/discordant_lines.txt | cut -f1 | sort -u | wc -l) reads differ" || echo yes)" | tee -a $D/summary.txt
awk '{for(i=14;i<=NF;i++) if($i~/^AS:i:/) print $1, $2, $7":"$9"-"$10, $i}' $D/discordant_lines.txt | tee -a $D/summary.txt
for rep in 1 2 3; do for c in stock heuristic certified; do
  case $c in stock) E="MM2_GEO=0";; heuristic) E="MM2_GEO=1 MM2_GEO_FAST=1";; certified) E="MM2_GEO=1";; esac
  /usr/bin/time -f "%e %M" -o $D/t env $E "$MM2GEO" -ax map-ont -t"$THREADS" grch38_ont.mmi $FQ > /dev/null 2>/dev/null
  echo -e "$c\trep$rep\t$(cat $D/t)" | tee -a $D/timing.tsv; done; done
for c in stock heuristic certified; do
  case $c in stock) E="MM2_GEO=0";; heuristic) E="MM2_GEO=1 MM2_GEO_FAST=1";; certified) E="MM2_GEO=1";; esac
  /usr/bin/time -f "%e %U %S" -o $D/t env $E MM2_GEO_TIMEDP=1 "$MM2GEO_INSTR" -ax map-ont -t"$THREADS" grch38_ont.mmi $FQ > /dev/null 2> $D/err
  read e u s < $D/t; echo -e "$c\twall=$e\tcpu=$(echo "$u + $s" | bc)\tdp_cpu=$(grep -o 'dp_cpu_s=[0-9.]*' $D/err | cut -d= -f2)" | tee -a $D/dpstage.tsv; done
