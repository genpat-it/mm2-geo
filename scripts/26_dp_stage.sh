#!/bin/bash
# DP-stage decomposition (Table 2 and Supplementary Table S21): thread-CPU time inside ksw2 (instrumented build,
# MM2_GEO_TIMEDP=1) and total CPU time for minimap2 / heuristic / certified; one run each at $THREADS threads.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/dpstage; mkdir -p $D
[ -s grch38_hifi.mmi ] || MM2_GEO=0 "$MM2GEO" -x map-hifi -d grch38_hifi.mmi "$GRCH38_REF" > /dev/null 2>&1
echo -e "data\tmode\twall_s\tcpu_s\tdp_cpu_s" > $D/dpstage.tsv
for s in "hifi_chr22 map-hifi chr22_named.fa hifi_chr22.fq" "ont_chr22 map-ont chr22_named.fa ont_r10_chr22.fq" "clr_chr22 map-pb chr22_named.fa clr_chr22.fq" "ont_chr14 map-ont ref_chr14.fa ont_r10_chr14.fq" "hifi_genomewide_subset map-hifi grch38_hifi.mmi $GW_READS"; do set -- $s
  for c in stock heuristic certified; do
    case $c in stock) E="MM2_GEO=0";; heuristic) E="MM2_GEO=1";; certified) E="MM2_GEO=1 MM2_GEO_CERT=1";; esac
    /usr/bin/time -f "%e %U %S" -o $D/t env $E MM2_GEO_TIMEDP=1 "$MM2GEO_INSTR" -ax $2 -t"$THREADS" $3 $4 > /dev/null 2> $D/err
    read e u s < $D/t; echo -e "$1\t$c\t$e\t$(echo "$u + $s" | bc)\t$(grep -o 'dp_cpu_s=[0-9.]*' $D/err | cut -d= -f2)" | tee -a $D/dpstage.tsv
  done
done
python3 - $D/dpstage.tsv <<'PY'
import sys, collections
r=collections.defaultdict(dict)
for l in list(open(sys.argv[1]))[1:]:
    d,m,w,c,dp=l.split('\t'); r[d][m]=(float(w),float(c),float(dp))
for d,v in r.items():
    w0,c0,p0=v['stock']
    print(f"{d}: DP share {100*p0/c0:.0f}%", *(f"| {m}: DP {p0/v[m][2]:.2f}x, non-DP {100*((v[m][1]-v[m][2])/(c0-p0)-1):+.1f}%, CPU {c0/v[m][1]:.2f}x, wall {w0/v[m][0]:.2f}x, predicted {c0/((c0-p0)+v[m][2]):.2f}x" for m in ('heuristic','certified')))
PY
