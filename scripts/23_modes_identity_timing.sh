#!/bin/bash
# Operating points (Table 2 end-to-end columns, Supplementary Table S22): output identity of the certified mode
# (PAF with CIGAR and SAM, sha256) and timing of minimap2 / heuristic / certified, same executable,
# 3 interleaved runs, SAM to /dev/null. chr22 HiFi/ONT/CLR and chr14 ONT.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/modes; mkdir -p $D
for s in "hifi_chr22 map-hifi chr22_named.fa hifi_chr22.fq" "ont_chr22 map-ont chr22_named.fa ont_r10_chr22.fq" "clr_chr22 map-pb chr22_named.fa clr_chr22.fq" "ont_chr14 map-ont ref_chr14.fa ont_r10_chr14.fq"; do set -- $s
  for fmt in c a; do
    h0=$(MM2_GEO=0 "$MM2GEO" -${fmt}x $2 -t"$THREADS" $3 $4 2>/dev/null | grep -v '^@PG' | sha256sum | cut -c1-16)
    h1=$(MM2_GEO=1 MM2_GEO_CERT=1 "$MM2GEO" -${fmt}x $2 -t"$THREADS" $3 $4 2>/dev/null | grep -v '^@PG' | sha256sum | cut -c1-16)
    echo -e "$1\t-${fmt}x\tcertified identical to minimap2: $([ $h0 = $h1 ] && echo yes || echo NO)" | tee -a $D/identity.tsv
  done
  for rep in 1 2 3; do for c in stock heuristic certified; do
    case $c in stock) E="MM2_GEO=0";; heuristic) E="MM2_GEO=1";; certified) E="MM2_GEO=1 MM2_GEO_CERT=1";; esac
    /usr/bin/time -f "%e %M" -o $D/t env $E "$MM2GEO" -ax $2 -t"$THREADS" $3 $4 > /dev/null 2>/dev/null
    echo -e "$1\t$c\trep$rep\t$(cat $D/t)" | tee -a $D/timing.tsv
  done; done
done
python3 - "$D/timing.tsv" <<'PY'
import sys, collections, statistics as st
t=collections.defaultdict(list)
for l in open(sys.argv[1]):
    d,c,r,x=l.rstrip('\n').split('\t'); t[(d,c)].append(float(x.split()[0]))
for d in dict.fromkeys(k[0] for k in t):
    s=st.median(t[(d,'stock')]); print(d, ' '.join(f"{c} {st.median(t[(d,c)]):.1f}s ({s/st.median(t[(d,c)]):.2f}x)" for c in ('stock','heuristic','certified')))
PY
