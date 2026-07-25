#!/bin/bash
# Alignment concordance mm2 vs geo (chr22, GRCh38): primary reads only.
# Reports %same-POS, %same-MAPQ, %same-CIGAR, and NM (edit-distance) / AS (score) deltas.
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"; REF=chr22_named.fa
SAM="$SAMTOOLS"
R=hifi_bench/results_CONCORDANCE.txt; : > $R
D=conc_scratch; mkdir -p $D
extract(){ # bam-stream -> name POS MAPQ CIGAR NM AS   (primary only)
  $SAM view -F 0x900 - 2>/dev/null | awk 'BEGIN{OFS="\t"}{nm="NA";as="NA";
    for(i=12;i<=NF;i++){if($i~/^NM:i:/){nm=substr($i,6)} if($i~/^AS:i:/){as=substr($i,6)}}
    print $1,$4,$5,$6,nm,as}' | sort -k1,1
}
run(){ local tech=$1 preset=$2 reads=$3
  echo ">> $tech concordance (mm2 vs geo)" | tee -a $R
  MM2_GEO=0 $MG -ax $preset -t32 $REF $reads 2>/dev/null | extract > $D/$tech.mm2
  MM2_GEO=1 MM2_GEO_MARGIN=20 $MG -ax $preset -t32 $REF $reads 2>/dev/null | extract > $D/$tech.geo
  join -t$'\t' $D/$tech.mm2 $D/$tech.geo | awk -F'\t' '
    {n++;
     if($2==$7)pos++; if($3==$8)mapq++; if($4==$9)cig++;
     if($5!="NA"&&$10!="NA"){dnm=$10-$5; anm+=(dnm<0?-dnm:dnm); if(dnm>0)nmworse++; else if(dnm<0)nmbetter++; nmn++}
     if($6!="NA"&&$11!="NA"){das=$11-$6; aas+=(das<0?-das:das); asn++}
    }
    END{printf "   shared primaries: %d\n   same POS: %.2f%%  same MAPQ: %.2f%%  identical CIGAR: %.2f%%\n",n,100*pos/n,100*mapq/n,100*cig/n;
        if(nmn)printf "   mean |dNM|: %.3f  (geo worse: %.2f%%, better: %.2f%%, equal: %.2f%%)\n",anm/nmn,100*nmworse/nmn,100*nmbetter/nmn,100*(nmn-nmworse-nmbetter)/nmn;
        if(asn)printf "   mean |dAS|: %.3f\n",aas/asn}' | tee -a $R
}
run HiFi map-hifi hifi_chr22.fq
run ONT  map-ont  ont_r10_chr22.fq
rm -rf $D; echo "### DONE-CONC" | tee -a $R
