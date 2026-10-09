#!/bin/bash
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"; SAM="$SAMTOOLS"
REF="$GRCH38_REF"
RD="$GW_READS"
R=hifi_bench/results_WG_CONCORD.txt; : > $R
export LC_ALL=C
# read stats
awk 'NR%4==2{n++; b+=length($0); print length($0)}' $RD | sort -rn > $TMPDIR/wglens.txt
awk -v tot=$(awk '{s+=$1}END{print s}' $TMPDIR/wglens.txt) 'BEGIN{h=tot/2}{c+=$1; if(c>=h && !done){print "N50",$1; done=1}} END{}' $TMPDIR/wglens.txt >> $R
awk 'END{}' ; nb=$(awk '{s+=$1}END{print s}' $TMPDIR/wglens.txt); nr=$(wc -l <$TMPDIR/wglens.txt)
awk -v nr=$nr -v nb=$nb 'BEGIN{printf "reads %d  bases %d (%.2f Gb)  mean_len %d  approx_cov %.2fx (vs 3.1Gb)\n",nr,nb,nb/1e9,nb/nr,nb/3.1e9}' | tee -a $R
D=wgc; mkdir -p $D
# primary: name, RNAME, strand, POS, MAPQ, CIGAR
ex(){ $SAM view -F 0x904 - 2>/dev/null | awk 'BEGIN{OFS="\t"}{s=and($2,16)?"-":"+"; print $1,$3,s,$4,$5,$6}' | LC_ALL=C sort -t$'\t' -k1,1; }
MM2_GEO=0 $MG -ax map-hifi -t32 $REF $RD 2>/dev/null | ex > $D/mm2
MM2_GEO=1 $MG -ax map-hifi -t32 $REF $RD 2>/dev/null | ex > $D/geo
nmm2=$(wc -l <$D/mm2); ngeo=$(wc -l <$D/geo)
join -t$'\t' $D/mm2 $D/geo > $D/j
awk -v a=$nmm2 -v b=$ngeo -F'\t' '
{sh++;
 if($2==$7)sc++; if($3==$8)sst++; if($4==$9)sp++; if($5==$10)mq++; if($6==$11)cg++;
 if($2==$7 && $3==$8 && $4==$9)place++}
END{printf "primary shared %d (mm2 %d, geo %d)\n",sh,a,b
 printf "same contig %.2f%%  same strand %.2f%%  same POS %.2f%%  same placement(contig+strand+POS) %.2f%%\n",100*sc/sh,100*sst/sh,100*sp/sh,100*place/sh
 printf "same MAPQ %.2f%%  CIGARid %.2f%%\n",100*mq/sh,100*cg/sh
 printf "diff contig %d  mm2-only %d  geo-only %d\n",sh-sc,a-sh,b-sh}' $D/j | tee -a $R
rm -rf $D; echo "### DONE-WGCONCORD" | tee -a $R
