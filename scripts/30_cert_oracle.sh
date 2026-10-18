#!/bin/bash
# Per-gap oracle for the certified mode (Supplementary, "Proof of the score certificate"): every gap the
# certificate accepts below the full band (and not z-dropped) is re-aligned score-only at the full band with
# z-drop off; a violation is a full-band score above the accepted score. Positive control: the same oracle on
# the heuristic mode (MM2_GEO_CERT_ORACLE=2), which is known to accept lower-scoring paths.
# Needs MM2GEO_ORACLE (revision-2026 + cert_oracle.patch; the hook does not change the output) and the
# stress sets of scripts 13 and 14 and the E. coli set of script 16.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/cert_oracle; mkdir -p $D
: > $D/summary.tsv
run(){ # name preset ref reads oracle_mode cert
  env MM2_GEO=1 MM2_GEO_CERT=$6 MM2_GEO_CERT_ORACLE=$5 "$MM2GEO_ORACLE" -cx $2 -t"$THREADS" $3 $4 2> >(grep '^ORC' | gzip > $D/$1.orc.gz) > /dev/null
  wait; sleep 2
  zcat $D/$1.orc.gz | awk -v n=$1 '{g++; if($10){z++; next} if($9==$8) eq++; if($9>$8){v++; if($9-$8>mx) mx=$9-$8}}
    END{printf "%s\tgaps_checked=%d\tviolations=%d\tmax_excess=%d\tequal_score=%.4f%%\tfull_band_skipped=%d\n", n, g, v+0, mx+0, 100*eq/(g-z), z+0}' | tee -a $D/summary.tsv
}
run control_heuristic_hifi_chr22 map-hifi chr22_named.fa hifi_chr22.fq 2 0
run hifi_chr22    map-hifi chr22_named.fa hifi_chr22.fq     1 1
run ont_r10_chr22 map-ont  chr22_named.fa ont_r10_chr22.fq  1 1
run clr_chr22     map-pb   chr22_named.fa clr_chr22.fq      1 1
run dense_indels  map-hifi chr22_named.fa adv_chr22/adv.fq  1 1
run large_sv      map-hifi chr22_named.fa advsv_chr22/sv.fq 1 1
run ecoli_ont     map-ont  ecoli_vc/ref_mut.fa ecoli_real2/SRR9900640.fastq 1 1
