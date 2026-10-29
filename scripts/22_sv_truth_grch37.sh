#!/bin/bash
# SV calling against GIAB HG002 Tier1 v0.6 (GRCh37 chr22) for minimap2, heuristic and certified modes, HiFi/ONT/CLR
# (extends 07_sv_calling.sh with the certified mode and CLR; Supplementary Table S17, bottom block).
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; S=sv; D=results/sv_truth; mkdir -p $D
: > $D/summary.txt   # fresh output on every run (no mixing with a previous run)
REF=$S/grch37_chr22.fa; CTG=$(head -1 $REF | sed 's/>//;s/ .*//')
[ -s $S/truth_$CTG.vcf.gz ] || { "$BCFTOOLS" view $S/HG002_SVs_Tier1_v0.6.vcf.gz -r $CTG -Oz -o $S/truth_$CTG.vcf.gz; "$BCFTOOLS" index -tf $S/truth_$CTG.vcf.gz; }
[ -s $S/tier1_$CTG.bed ] || awk -v c=$CTG '$1==c' $S/HG002_SVs_Tier1_v0.6.bed > $S/tier1_$CTG.bed
for s in "HiFi map-hifi hifi_chr22.fq" "ONT map-ont ont_r10_chr22.fq" "CLR map-pb clr_chr22.fq"; do set -- $s
  for c in stock heuristic certified; do
    case $c in stock) E="MM2_GEO=0";; heuristic) E="MM2_GEO=1 MM2_GEO_FAST=1";; certified) E="MM2_GEO=1";; esac
    b=$D/${1}_$c.bam
    [ -s $b.bai ] || { env $E "$MM2GEO" -ax $2 -t"$THREADS" $REF $3 2>/dev/null | "$SAMTOOLS" sort -@8 -m2G -o $b - 2>/dev/null; "$SAMTOOLS" index $b; }
    "$SNIFFLES" --input $b --vcf $D/${1}_$c.vcf.gz --reference $REF --threads 16 > $D/snf_${1}_$c.log 2>&1
    rm -rf $D/tv_${1}_$c; "$TRUVARI" bench -b $S/truth_$CTG.vcf.gz -c $D/${1}_$c.vcf.gz --includebed $S/tier1_$CTG.bed -o $D/tv_${1}_$c > /dev/null 2>&1
    python3 -c "import json;d=json.load(open('$D/tv_${1}_$c/summary.json'));print('$1 $c P=%.3f R=%.3f F1=%.3f TP=%d FP=%d FN=%d'%(d['precision'],d['recall'],d['f1'],d['TP-comp'],d['FP'],d['FN']))" | tee -a $D/summary.txt
    echo "   615 bp VNTR insertion HG2_PB_SVrefine2Falcon1plusDovetail_7697 recovered: $(zcat $D/tv_${1}_$c/tp-base.vcf.gz | grep -c HG2_PB_SVrefine2Falcon1plusDovetail_7697)" | tee -a $D/summary.txt
    # read-level view of the 615 bp insertion (GRCh37 chr22:47,566,242): reads with a primary or supplementary
    # alignment within 300 bp of the site, reads whose alignment over the site carries one insertion >= 300 bp,
    # and the Sniffles2 call there
    echo "   reads within 300 bp of the site $("$SAMTOOLS" view -F 0x904 $b $CTG:47565942-47566542 | cut -f1 | sort -u | wc -l), with a single insertion >=300 bp $("$SAMTOOLS" view -F 0x904 $b $CTG:47566242-47566242 | awk '{c=$6; while (match(c,/^[0-9]+[MIDNSHP=X]/)) { if (substr(c,RLENGTH,1)=="I" && substr(c,RSTART,RLENGTH-1)+0>=300) {print $1; break}; c=substr(c,RLENGTH+1) }}' | sort -u | wc -l)" | tee -a $D/summary.txt
    zcat $D/${1}_$c.vcf.gz | awk -v s=47566242 '!/^#/ && $8~/SVTYPE=INS/ && $2>s-500 && $2<s+500 {
        match($8,/SVLEN=-?[0-9]+/); l=substr($8,RSTART+6,RLENGTH-6); match($8,/SUPPORT=[0-9]+/); r=substr($8,RSTART+8,RLENGTH-8)
        printf "   Sniffles2 call at the site: INS %s bp, %s reads\n", l, r; f=1} END{if(!f) print "   Sniffles2 call at the site: none"}' | tee -a $D/summary.txt
  done
done
