#!/bin/bash
# Whole-genome SPEED (alignment-only): real genome-wide HG002 HiFi -> full GRCh38, mm2 vs geo.
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"; SAM="$SAMTOOLS"
REF="$GRCH38_REF"
RD="$GW_READS"
export LC_ALL=C; R=hifi_bench/results_WHOLEGENOME.txt; : > $R; D=wg; mkdir -p $D
echo "### WHOLE-GENOME speed (real genome-wide HG002 HiFi $(($(wc -l <$RD)/4)) reads -> GRCh38 full), mm2 vs geo, 32t, one untimed warm-up run, then median of 3 interleaved runs" | tee -a $R
[ -s $D/grch38.mmi ] || { echo "building index..."|tee -a $R; MM2_GEO=0 $MG -x map-hifi -d $D/grch38.mmi $REF >/dev/null 2>&1; }
# warm-up: loads the 7 GB index and the reads into the page cache, so that no timed run pays for cold storage
MM2_GEO=0 $MG -ax map-hifi -t32 $D/grch38.mmi $RD >/dev/null 2>/dev/null
: > $D/times
for rep in 1 2 3; do for cfg in mm2 geo; do
  gf=0; [ "$cfg" = geo ] && gf=1
  /usr/bin/time -f "%e %M" bash -c "MM2_GEO=$gf MM2_GEO_MARGIN=20 $MG -ax map-hifi -t32 $D/grch38.mmi $RD >/dev/null 2>/dev/null" 2>$D/t
  echo "$cfg $(cat $D/t)" >> $D/times
done; done
for cfg in mm2 geo; do
  awk -v c=$cfg '$1==c{print $2, $3}' $D/times | sort -n | awk -v c=$cfg 'NR==2{printf "   %-4s time %ss  RSS %.1fGB (median of 3)\n",c,$1,$2/1048576}' | tee -a $R
done
ex(){ $SAM view -F 0x900 - 2>/dev/null | awk 'BEGIN{OFS="\t"}{nm="NA";for(i=12;i<=NF;i++)if($i~/^NM:i:/)nm=substr($i,6);print $1,$4,$5,$6,nm}'|sort -k1,1; }
MM2_GEO=0 $MG -ax map-hifi -t32 $D/grch38.mmi $RD 2>/dev/null | tee >($SAM view -c -F 0x900 - 2>/dev/null >$D/nmap_mm2) | ex >$D/mm2
MM2_GEO=1 MM2_GEO_MARGIN=20 $MG -ax map-hifi -t32 $D/grch38.mmi $RD 2>/dev/null | ex >$D/geo
echo "   primary alignments: mm2 $(cat $D/nmap_mm2 2>/dev/null), geo $(wc -l <$D/geo)" | tee -a $R
join -t$'\t' $D/mm2 $D/geo | awk -F'\t' '{n++;if($2==$6)p++;if($3==$7)mq++;if($4==$8)cg++;if($5!="NA"&&$9!="NA"&&$5==$9)nm++}
  END{printf "   concordance (%d shared): samePOS %.2f%% sameMAPQ %.2f%% CIGARid %.2f%% NMid %.2f%%\n",n,100*p/n,100*mq/n,100*cg/n,100*nm/n}' | tee -a $R
rm -rf $D/mm2 $D/geo; echo "### DONE-WG" | tee -a $R
