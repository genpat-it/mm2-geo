#!/bin/bash
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"; SAM="$SAMTOOLS"
export LC_ALL=C TMPDIR="$WORKDIR"/prtmp
mkdir -p $TMPDIR; R=hifi_bench/results_PRESENCE.txt; : > $R
names(){ $SAM view -F 0x904 - 2>/dev/null | cut -f1 | LC_ALL=C sort -T $TMPDIR -u; }
row(){ local lab=$1 pre=$2 rd=$3 ref=$4
  MM2_GEO=0 $MG -ax $pre -t32 $ref $rd 2>/dev/null | names > $TMPDIR/m
  MM2_GEO=1 $MG -ax $pre -t32 $ref $rd 2>/dev/null | names > $TMPDIR/g
  local nm=$(wc -l <$TMPDIR/m) ng=$(wc -l <$TMPDIR/g)
  local sh=$(comm -12 $TMPDIR/m $TMPDIR/g | wc -l)
  local mo=$(comm -23 $TMPDIR/m $TMPDIR/g | wc -l)
  local go=$(comm -13 $TMPDIR/m $TMPDIR/g | wc -l)
  printf "%-12s mm2_primary %d  geo_primary %d  shared %d  mm2_only %d  geo_only %d\n" "$lab" $nm $ng $sh $mo $go | tee -a $R
}
row "chr22-HiFi" map-hifi hifi_chr22.fq chr22_named.fa
row "chr22-ONT"  map-ont  ont_r10_chr22.fq chr22_named.fa
row "chr22-CLR"  map-pb   clr_chr22.fq chr22_named.fa
row "chr14-CLR"  map-pb   clr_chr14.fq ref_chr14.fa
row "ecoli-ONT"  map-ont  ecoli_real2/SRR9900640.fastq "$ECOLI_REF"
rm -rf $TMPDIR; echo "### DONE-PRESENCE" | tee -a $R
