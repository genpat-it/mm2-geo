#!/bin/bash
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"; SAM="$SAMTOOLS"
AR=chr22_named.fa; R=hifi_bench/results_ABLCIGAR_FULL.txt; : > $R
export LC_ALL=C
# same extraction as run_concordance.sh: primary only (-F 0x900), name<TAB>CIGAR
ex(){ $SAM view -F 0x900 - 2>/dev/null | awk 'BEGIN{OFS="\t"}{print $1,$6}' | sort -k1,1; }
D=abcf; mkdir -p $D
for TE in "hifi map-hifi hifi_chr22.fq" "ont map-ont ont_r10_chr22.fq"; do
  set -- $TE; tech=$1; pre=$2; rd=$3
  MM2_GEO=0 $MG -ax $pre -t32 $AR $rd 2>/dev/null | ex > $D/base
  for cfg in "fixed-narrow|MM2_GEO=1 MM2_GEO_FIXED=1 MM2_GEO_NODOUBLE=1" \
             "fixed-doubling|MM2_GEO=1 MM2_GEO_FIXED=1" \
             "geometry-only|MM2_GEO=1 MM2_GEO_NODOUBLE=1" \
             "mm2-geo|MM2_GEO=1"; do
    name="${cfg%%|*}"; env="${cfg#*|}"
    eval "$env $MG -ax $pre -t32 $AR $rd 2>/dev/null" | ex > $D/cfg
    join -t$'\t' $D/base $D/cfg | awk -v t=$tech -v n="$name" -F'\t' '{c++; if($2==$3)m++}END{printf "%-5s %-16s shared %d CIGARid %.2f%%\n",t,n,c,100*m/c}' | tee -a $R
  done
done
rm -rf $D; echo "### DONE-ABLCIGAR" | tee -a $R
