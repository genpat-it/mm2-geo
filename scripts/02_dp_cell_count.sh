#!/bin/bash
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"
R=/tmp/mechdef.out; : > $R
for T in "HiFi map-hifi /tmp/c22.hifi.mmi hifi_chr22.fq" "ONT map-ont /tmp/c22.ont.mmi ont_r10_chr22.fq"; do
  set -- $T; tech=$1; pre=$2; idx=$3; rd=$4
  echo "==== $tech (full chr22) ====" | tee -a $R
  # mm2 cells
  mc=$(MM2_GEO=0 MM2_GEO_STATS=1 $MG -ax $pre -t32 $idx $rd 2>&1 >/dev/null | grep MM2_GEO_CELLS | grep -o 'dp_cells=[0-9]*' | cut -d= -f2)
  # geo cells + band stats (one run)
  MM2_GEO=1 MM2_GEO_STATS=1 $MG -ax $pre -t32 $idx $rd 2>/tmp/geo_$tech.err >/dev/null
  gc=$(grep MM2_GEO_CELLS /tmp/geo_$tech.err | grep -o 'dp_cells=[0-9]*' | cut -d= -f2)
  st=$(grep MM2_GEO_STATS /tmp/geo_$tech.err)
  awk -v t=$tech -v m=$mc -v g=$gc 'BEGIN{printf "%s  mm2_cells=%s  geo_cells=%s  rel=%.1f%%  reduction=%.2fx\n",t,m,g,(m?100*g/m:0),(g?m/g:0)}' | tee -a $R
  echo "  $st" | tee -a $R
done
echo "### MECHDEF-DONE" | tee -a $R
