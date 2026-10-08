#!/bin/bash
# Full-coverage genome-wide HG002 HiFi benchmark per partition (Supplementary Table S18): reads partitioned by
# chromosome of origin (primary records of the GIAB alignment), each partition aligned to the complete GRCh38
# (pre-built index), 3 interleaved runs. Heuristic series: submission-2026 executable; certified series:
# revision-2026 executable, with a PAF identity check per partition. Also the chr1-only setting.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/genomewide; mkdir -p $D tmp_wg
IDX=grch38_hifi.mmi; [ -s $IDX ] || MM2_GEO=0 "$MM2GEO" -x map-hifi -d $IDX "$GRCH38_REF" > /dev/null 2>&1
[ -s "$HG002_HIFI_BAM" ] || { echo "run: bash scripts/fetch_data.sh genomewide"; exit 1; }
"$SAMTOOLS" idxstats "$HG002_HIFI_BAM" | awk '$1!~/^chr([0-9]+|X|Y)$/ && $1!="*"{print $1}' > tmp_wg/other.lst
for p in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY other; do
  grep -q "^$p"$'\t' $D/perchr.tsv 2>/dev/null && continue
  fq=tmp_wg/$p.fq; reg=$p; [ $p = other ] && reg="$(cat tmp_wg/other.lst)"
  "$SAMTOOLS" view -@8 -b -F 0x900 "$HG002_HIFI_BAM" $reg | "$SAMTOOLS" fastq -@8 - > $fq 2>/dev/null
  nr=$(awk 'NR%4==1{n++}END{print n+0}' $fq)
  a=$(MM2_GEO=0 "$MM2GEO" -cx map-hifi -t"$THREADS" $IDX $fq 2>/dev/null | sort -S4G | sha256sum | cut -c1-16)
  b=$(MM2_GEO=1 MM2_GEO_CERT=1 "$MM2GEO" -cx map-hifi -t"$THREADS" $IDX $fq 2>/dev/null | sort -S4G | sha256sum | cut -c1-16)
  idn=$([ $a = $b ] && echo yes || echo no)
  for rep in 1 2 3; do
    for c in stock_sub heuristic stock certified; do
      case $c in stock_sub) X="$MM2GEO_SUB"; E="MM2_GEO=0";; heuristic) X="$MM2GEO_SUB"; E="MM2_GEO=1";; stock) X="$MM2GEO"; E="MM2_GEO=0";; certified) X="$MM2GEO"; E="MM2_GEO=1 MM2_GEO_CERT=1";; esac
      /usr/bin/time -f "%e %M" -o tmp_wg/t env $E "$X" -ax map-hifi -t"$THREADS" $IDX $fq > /dev/null 2>/dev/null
      echo -e "$p\t$nr\t$idn\t$c\t$rep\t$(cat tmp_wg/t)" >> $D/perchr.tsv
    done
  done
  [ $p = chr1 ] && { "$SAMTOOLS" faidx "$GRCH38_REF" CM000663.2 > tmp_wg/chr1.fa; "$MM2GEO" -x map-hifi -d tmp_wg/chr1.mmi tmp_wg/chr1.fa > /dev/null 2>&1
    for rep in 1 2 3; do for g in 0 1; do /usr/bin/time -f "%e %M" -o tmp_wg/t env MM2_GEO=$g "$MM2GEO_SUB" -ax map-hifi -t"$THREADS" tmp_wg/chr1.mmi $fq > /dev/null 2>/dev/null
      echo -e "chr1_only_reference\tgeo=$g\trep$rep\t$(cat tmp_wg/t)" >> $D/chr1_only.tsv; done; done; }
  rm -f $fq
done
python3 - $D/perchr.tsv <<'PY'
import sys, collections, statistics as st
t=collections.defaultdict(list); idn={}
for l in open(sys.argv[1]):
    p,n,i,c,r,x=l.rstrip('\n').split('\t'); t[(p,c)].append(float(x.split()[0])); idn[p]=i
tot=collections.Counter()
for p in dict.fromkeys(k[0] for k in t):
    m={c:st.median(t[(p,c)]) for c in ('stock_sub','heuristic','stock','certified')}; tot.update(m)
    print(f"{p}: heuristic {m['stock_sub']/m['heuristic']:.2f}x  certified {m['stock']/m['certified']:.2f}x  certified identical: {idn[p]}")
print(f"whole genome: heuristic {tot['stock_sub']/tot['heuristic']:.2f}x  certified {tot['stock']/tot['certified']:.2f}x")
PY
