#!/bin/bash
# Certified-mode counterparts of the experiments first run with the heuristic mode only. Same executable (tag
# revision-2026) for minimap2 (MM2_GEO=0), the certified mode (MM2_GEO=1) and the heuristic mode
# (MM2_GEO=1 MM2_GEO_FAST=1); configurations interleaved within each replicate. Output: results/certified/.
#   chr1only   full-coverage chr1 reads against a chr1-only index (Table S18 note), 32 threads, median of 3
#   speed      Table S20 settings: 1/4/8/32 threads, SAM to /dev/null or disk, chr1-only and GRCh38 index,
#              -O3 rebuild and bioconda minimap2 (genome-wide subset reads; 30k reads at 1 thread)
#   threads    thread scaling 1-64 (HiFi chr22, 100k reads; single run at 1-2 threads, median of 3 at >=4),
#              with mm2-fast and mm2-fast+geo in both modes
#   margin     margin m = 5..80 (chr22 HiFi/ONT), time of both modes and identity of the certified output to minimap2
#   rband      minimap2 -r 20,20 vs default vs both modes (chr22 HiFi/ONT), time and certified identity
#   gwsub      genome-wide subset (0.57x) against GRCh38 (pre-built index, one warm-up run), median of 3
#   cells      logical DP cells (MM2_GEO_STATS), chr22 HiFi/ONT/CLR
#   ecoli      certified output identical to minimap2 on the E. coli isolate (both references)
# Usage: bash 37_certified_complements.sh [section ...]   (default: all sections)
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/certified; mkdir -p $D
SECTIONS="${*:-chr1only speed threads margin rband gwsub cells ecoli}"
E_mm2="MM2_GEO=0"; E_cert="MM2_GEO=1"; E_heur="MM2_GEO=1 MM2_GEO_FAST=1"
tm(){ # tsv label env executable args...   (time + RSS of one run, SAM to /dev/null)
  local f=$1 l=$2 e=$3 x=$4; shift 4
  /usr/bin/time -f "%e %M" -o $D/t env $e "$x" "$@" > /dev/null 2>/dev/null; echo -e "$l\t$(tr ' ' '\t' < $D/t)" >> $f; }
same(){ # preset ref reads -> yes/NO (certified SAM equals minimap2's, @PG excluded)
  local a b; a=$(MM2_GEO=0 "$MM2GEO" -ax $1 -t32 $2 $3 2>/dev/null | grep -v '^@PG' | md5sum | cut -c1-32)
  b=$(MM2_GEO=1 "$MM2GEO" -ax $1 -t32 $2 $3 2>/dev/null | grep -v '^@PG' | md5sum | cut -c1-32); [ "$a" = "$b" ] && echo yes || echo NO; }
for sec in $SECTIONS; do case $sec in
chr1only)
  f=$D/chr1only.tsv; : > $f; mkdir -p tmp_wg
  [ -s tmp_wg/chr1.mmi ] || { "$SAMTOOLS" faidx "$GRCH38_REF" CM000663.2 > tmp_wg/chr1.fa; MM2_GEO=0 "$MM2GEO" -x map-hifi -d tmp_wg/chr1.mmi tmp_wg/chr1.fa > /dev/null 2>&1; }
  "$SAMTOOLS" view -@8 -b -F 0x900 "$HG002_HIFI_BAM" chr1 | "$SAMTOOLS" fastq -@8 - > tmp_wg/chr1.fq 2>/dev/null
  for rep in 1 2 3; do for c in mm2 cert heur; do v=E_$c; tm $f "$c\t$rep" "${!v}" "$MM2GEO" -ax map-hifi -t32 tmp_wg/chr1.mmi tmp_wg/chr1.fq; done; done
  echo -e "identity\t$(same map-hifi tmp_wg/chr1.mmi tmp_wg/chr1.fq)" >> $f; rm -f tmp_wg/chr1.fq ;;
speed)
  f=$D/speed.tsv; : > $f; S=$D/speed_tmp; mkdir -p $S
  [ -s $S/chr1.mmi ] || { "$SAMTOOLS" faidx "$GRCH38_REF" CM000663.2 > $S/chr1.fa; MM2_GEO=0 "$MM2GEO" -x map-hifi -d $S/chr1.mmi $S/chr1.fa > /dev/null 2>&1; }
  [ -s grch38_hifi.mmi ] || MM2_GEO=0 "$MM2GEO" -x map-hifi -d grch38_hifi.mmi "$GRCH38_REF" > /dev/null 2>&1
  head -n 120000 "$GW_READS" > $S/reads_30k.fq
  for t in 1 4 8 32; do rd="$GW_READS"; [ $t = 1 ] && rd=$S/reads_30k.fq
    for out in devnull disk; do
      for rep in 1 2 3; do for c in mm2 cert; do v=E_$c
        if [ $out = disk ]; then /usr/bin/time -f "%e %M" -o $D/t env ${!v} "$MM2GEO" -ax map-hifi -t$t $S/chr1.mmi "$rd" > $S/out.sam 2>/dev/null
             echo -e "t$t\t$out\t$c\t$rep\t$(tr ' ' '\t' < $D/t)" >> $f; rm -f $S/out.sam
        else tm $f "t$t\t$out\t$c\t$rep" "${!v}" "$MM2GEO" -ax map-hifi -t$t $S/chr1.mmi "$rd"; fi
      done; done; done; done
  for rep in 1 2 3; do for c in mm2 cert; do v=E_$c; tm $f "grch38\tdevnull\t$c\t$rep" "${!v}" "$MM2GEO" -ax map-hifi -t32 grch38_hifi.mmi "$GW_READS"; done; done
  for t in 4 32; do for rep in 1 2 3; do
    tm $f "bioconda_t$t\tdevnull\tmm2\t$rep" "MM2_GEO=0" "$BIOCONDA_MM2" -ax map-hifi -t$t $S/chr1.mmi "$GW_READS"
    for c in mm2 cert; do v=E_$c; tm $f "O3_t$t\tdevnull\t$c\t$rep" "${!v}" "$MM2GEO_O3" -ax map-hifi -t$t $S/chr1.mmi "$GW_READS"; done; done; done
  rm -rf $S ;;
threads)
  f=$D/threads.tsv; : > $f; S=$D/ts_tmp; mkdir -p $S; head -400000 hifi_chr22.fq > $S/sub100k.fq
  for thr in 1 2 4 8 16 32 64; do n=3; [ $thr -le 2 ] && n=1
    for rep in $(seq 1 $n); do
      for c in mm2 heur cert; do v=E_$c; tm $f "$thr\t$c\t$rep" "${!v}" "$MM2GEO" -ax map-hifi -t$thr chr22_named.fa $S/sub100k.fq; done
      tm $f "$thr\tfast\t$rep" "" "$MM2FAST" -ax map-hifi -t$thr chr22_named.fa $S/sub100k.fq
      tm $f "$thr\tfastgeo_heur\t$rep" "$E_heur" "$MM2FASTGEO" -ax map-hifi -t$thr chr22_named.fa $S/sub100k.fq
      tm $f "$thr\tfastgeo_cert\t$rep" "$E_cert" "$MM2FASTGEO" -ax map-hifi -t$thr chr22_named.fa $S/sub100k.fq
    done; done
  rm -rf $S ;;
margin)
  f=$D/margin.tsv; : > $f
  for s in "HiFi map-hifi hifi_chr22.fq" "ONT map-ont ont_r10_chr22.fq"; do set -- $s
    ref=$(MM2_GEO=0 "$MM2GEO" -ax $2 -t32 chr22_named.fa $3 2>/dev/null | grep -v '^@PG' | md5sum | cut -c1-32)
    for rep in 1 2 3; do tm $f "$1\tmm2\t-\t$rep" "$E_mm2" "$MM2GEO" -ax $2 -t32 chr22_named.fa $3
      for m in 5 10 20 40 80; do for c in cert heur; do v=E_$c
        tm $f "$1\t$c\t$m\t$rep" "${!v} MM2_GEO_MARGIN=$m" "$MM2GEO" -ax $2 -t32 chr22_named.fa $3; done; done; done
    for m in 5 10 20 40 80; do
      x=$(MM2_GEO=1 MM2_GEO_MARGIN=$m "$MM2GEO" -ax $2 -t32 chr22_named.fa $3 2>/dev/null | grep -v '^@PG' | md5sum | cut -c1-32)
      echo -e "identity\t$1\t$m\t$([ "$x" = "$ref" ] && echo yes || echo NO)" >> $f; done
  done ;;
rband)
  f=$D/rband.tsv; : > $f
  for s in "HiFi map-hifi hifi_chr22.fq 20" "ONT map-ont ont_r10_chr22.fq 50"; do set -- $s
    for rep in 1 2 3; do
      tm $f "$1\tdefault\t$rep" "$E_mm2" "$MM2GEO" -ax $2 -t32 chr22_named.fa $3
      tm $f "$1\tr20\t$rep" "$E_mm2" "$MM2GEO" -ax $2 -r 20,20 -t32 chr22_named.fa $3
      tm $f "$1\theur\t$rep" "$E_heur MM2_GEO_MARGIN=$4" "$MM2GEO" -ax $2 -t32 chr22_named.fa $3
      tm $f "$1\tcert\t$rep" "$E_cert" "$MM2GEO" -ax $2 -t32 chr22_named.fa $3
    done; echo -e "identity\t$1\t$(same $2 chr22_named.fa $3)" >> $f; done ;;
gwsub)
  f=$D/gwsub.tsv; : > $f
  [ -s grch38_hifi.mmi ] || MM2_GEO=0 "$MM2GEO" -x map-hifi -d grch38_hifi.mmi "$GRCH38_REF" > /dev/null 2>&1
  MM2_GEO=0 "$MM2GEO" -ax map-hifi -t32 grch38_hifi.mmi "$GW_READS" > /dev/null 2>&1   # warm-up
  for rep in 1 2 3; do for c in mm2 heur cert; do v=E_$c; tm $f "$c\t$rep" "${!v}" "$MM2GEO" -ax map-hifi -t32 grch38_hifi.mmi "$GW_READS"; done; done
  echo -e "identity\t$(same map-hifi grch38_hifi.mmi "$GW_READS")" >> $f ;;
cells)
  f=$D/cells.tsv; : > $f
  for s in "HiFi map-hifi hifi_chr22.fq" "ONT map-ont ont_r10_chr22.fq" "CLR map-pb clr_chr22.fq"; do set -- $s
    for c in mm2 heur cert; do v=E_$c
      env ${!v} MM2_GEO_STATS=1 "$MM2GEO" -ax $2 -t32 chr22_named.fa $3 2>$D/st.err > /dev/null
      echo -e "$1\t$c\t$(grep -h 'MM2_GEO_CELLS\|MM2_GEO_STATS' $D/st.err | tr '\n' ' ')" >> $f; done; done; rm -f $D/st.err ;;
ecoli)
  f=$D/ecoli.tsv; : > $f
  echo -e "spikein_reference\t$(same map-ont ecoli_vc/ref_mut.fa ecoli_real2/SRR9900640.fastq)" >> $f
  echo -e "MG1655\t$(same map-ont "$ECOLI_REF" ecoli_real2/SRR9900640.fastq)" >> $f ;;
esac; echo "$sec done" >> $D/progress.txt; done
echo DONE >> $D/progress.txt
