#!/bin/bash
# UNIFIED CARTESIAN BENCHMARK: {dataset x tech x config} x {median-3 time, median-3 peak-RSS, F1}.
# Datasets: human chr22, human chr14 (GIAB real); E. coli (planted variants). Tech: HiFi/ONT/CLR.
# Configs: mm2 (stock), geo, mm2-fast, mm2-fast+geo. F1 via clair3 + rtg vcfeval.
set -uo pipefail
source "$(dirname "$0")/config.sh"
cd "$WORKDIR"
MG="$MM2GEO"; MF="$MM2FAST"; MFG="$MM2FASTGEO"
ECOLIREF="$ECOLI_REF"
C3="$CLAIR3"
M_HIFI="$MODEL_HIFI"
M_ONT="$MODEL_ONT"
M_EC_ONT="$MODEL_BACT"
source "$CONDA_PROFILE"
export LC_ALL=C; q(){ "$@" 2>/dev/null; }
S=cart_scratch; mkdir -p $S
TRES=hifi_bench/results_CART_timemem.txt; FRES=hifi_bench/results_CART_f1.txt; : > $TRES; : > $FRES
med(){ printf '%s\n' "$@"|sort -n|awk '{a[NR]=$1}END{print a[int((NR+1)/2)]}'; }

alncmd(){ case $1 in
  mm2)     echo "MM2_GEO=0 $MG -ax $2 -t32 $3 $4";;
  geo)     echo "MM2_GEO=1 MM2_GEO_MARGIN=20 $MG -ax $2 -t32 $3 $4";;
  fast)    echo "$MF -ax $2 -t32 $3 $4";;
  fastgeo) echo "MM2_GEO=1 MM2_GEO_MARGIN=20 $MFG -ax $2 -t32 $3 $4";;
 esac; }
timrss(){ local cmd; cmd="$(alncmd "$1" "$2" "$3" "$4")"; local t=() r=()
  for i in 1 2 3; do /usr/bin/time -f "%e %M" bash -c "$cmd >/dev/null 2>/dev/null" 2>$S/tr.t
    read e m < $S/tr.t; t+=("$e"); r+=("$m"); done
  awk -v tt=$(med "${t[@]}") -v rr=$(med "${r[@]}") 'BEGIN{printf "%5.0fs %5.1fGB",tt,rr/1048576}'; }

# datasets: name|alignref|clairref|ctg|truth|bed|sdf|reheaderTag
DS=(
 "chr22|chr22_named.fa|chr22_cm.fa|CM000684.2|truth_chr22.vcf.gz|conf_chr22.bed|ref_chr22_sdf|chr22"
 "chr14|ref_chr14.fa|ref_chr14.fa|CM000676.2|truth_chr14.vcf.gz|conf_chr14.bed|ref_chr14_sdf|NONE"
 "ecoli|$ECOLIREF|$ECOLIREF|NC_000913.3|ecoli_bench/truth_ecoli.vcf.gz|NONE|ecoli_bench/ecoli_sdf|NONE"
)
# tech per dataset: dsname:techname:preset:reads:model:platform  (model NONE=>no F1)
TECH=(
 "chr22:HiFi:map-hifi:hifi_chr22.fq:$M_HIFI:hifi"
 "chr22:ONT:map-ont:ont_r10_chr22.fq:$M_ONT:ont"
 "chr14:HiFi:map-hifi:hifi_chr14.fq:$M_HIFI:hifi"
 "chr14:ONT:map-ont:ont_r10_chr14.fq:$M_ONT:ont"
 "ecoli:HiFi:map-hifi:ecoli_bench/hifi.fq:$M_HIFI:hifi"
 "ecoli:ONT:map-ont:ecoli_bench/ont.fq:$M_EC_ONT:ont"
 "ecoli:CLR:map-pb:ecoli_bench/clr.fq:NONE:none"
)

echo "### PART A — TIME + PEAK-MEM (median of 3), all 4 configs" | tee -a $TRES
printf "%-6s %-5s | %-14s %-14s %-14s %-14s\n" "data" "tech" "mm2" "geo" "mm2-fast" "mm2-fast+geo" | tee -a $TRES
for t in "${TECH[@]}"; do
  IFS=: read dsn tech preset reads model plat <<< "$t"
  aref=""; for d in "${DS[@]}"; do IFS='|' read n ar cr ct tr bd sf rh <<< "$d"; [ "$n" = "$dsn" ] && aref=$ar; done
  a=$(timrss mm2 $preset $aref $reads); b=$(timrss geo $preset $aref $reads)
  c=$(timrss fast $preset $aref $reads); e=$(timrss fastgeo $preset $aref $reads)
  printf "%-6s %-5s | %-14s %-14s %-14s %-14s\n" "$dsn" "$tech" "$a" "$b" "$c" "$e" | tee -a $TRES
done
echo "### DONE-A" | tee -a $TRES

echo "### PART B — F1 vs truth (SNV/INDEL). mm2 & geo everywhere; +fast,fastgeo on chr22" | tee -a $FRES
tc(){ zcat "$1" 2>/dev/null | awk '!/^#/{if(length($4)==1&&length($5)==1)s++;else i++}END{printf "%d %d",s+0,i+0}'; }
dobam(){ # $1 cfg $2 preset $3 alignref $4 reads $5 reheaderFrom $6 reheaderTo $7 out.bam
  local cmd; cmd="$(alncmd "$1" "$2" "$3" "$4")"
  eval "$cmd" 2>/dev/null | q "$SAMTOOLS" sort -@8 -o $S/tmp.bam -
  if [ "$5" != "NONE" ]; then q "$SAMTOOLS" reheader -c "sed \"s/SN:$5/SN:$6/\"" $S/tmp.bam > $7; rm -f $S/tmp.bam; else mv $S/tmp.bam $7; fi
  "$SAMTOOLS" index -@8 $7 2>/dev/null; }
f1row(){ # label bam clairref ctg model plat truth bed sdf
  local lab=$1 bam=$2 cref=$3 ctg=$4 mdl=$5 plat=$6 tr=$7 bd=$8 sf=$9
  conda activate "$CONDA_ENV_CLAIR3"; rm -rf $S/c3
  local bedarg=""; [ "$bd" != "NONE" ] && bedarg="--bed_fn=$bd" || bedarg="--include_all_ctgs"
  $C3 --bam_fn=$bam --ref_fn=$cref --threads=32 --platform=$plat --model_path=$mdl --output=$S/c3 --ctg_name=$ctg $bedarg >$S/c3.log 2>&1
  conda deactivate; conda activate "$CONDA_ENV_RTG"; rm -rf $S/ev
  local sqarg=""; [ "$bd" = "NONE" ] && sqarg="--squash-ploidy"
  local ebed=""; [ "$bd" != "NONE" ] && ebed="-e $bd"
  rtg vcfeval -b $tr -c $S/c3/merge_output.vcf.gz -t $sf $ebed $sqarg -o $S/ev >/dev/null 2>&1 || { echo "  $lab FAIL"|tee -a $FRES; conda deactivate; return; }
  read ts ti < <(tc $S/ev/tp.vcf.gz); read fps fpi < <(tc $S/ev/fp.vcf.gz); read fns fni < <(tc $S/ev/fn.vcf.gz); read tbs tbi < <(tc $S/ev/tp-baseline.vcf.gz)
  conda deactivate
  awk -v ts=$ts -v ti=$ti -v fps=$fps -v fpi=$fpi -v fns=$fns -v fni=$fni -v tbs=$tbs -v tbi=$tbi -v l="$lab" \
   'function f(tp,fp,tb,fn, p,r){p=(tp+fp)?tp/(tp+fp):0;r=(tb+fn)?tb/(tb+fn):0;return sprintf("%.4f/%.4f/%.4f",p,r,(p+r)?2*p*r/(p+r):0)}
    BEGIN{printf "  %-22s SNV %-22s INDEL %-22s\n",l,f(ts,fps,tbs,fns),f(ti,fpi,tbi,fni)}' | tee -a $FRES; }
for t in "${TECH[@]}"; do
  IFS=: read dsn tech preset reads model plat <<< "$t"
  [ "$model" = "NONE" ] && { echo "  ${dsn} ${tech}: no caller model (CLR) — time/mem only" | tee -a $FRES; continue; }
  for d in "${DS[@]}"; do IFS='|' read n ar cr ct tr bd sf rh <<< "$d"; [ "$n" = "$dsn" ] && { AR=$ar;CR=$cr;CT=$ct;TR=$tr;BD=$bd;SF=$sf;RH=$rh; }; done
  cfgs="mm2 geo"; [ "$dsn" = "chr22" ] && cfgs="mm2 geo fast fastgeo"
  for cfg in $cfgs; do
    dobam $cfg $preset $AR $reads $RH $CT $S/b.bam
    f1row "${dsn} ${tech} ${cfg}" $S/b.bam $CR $CT $model $plat $TR $BD $SF
  done
done
rm -rf $S; echo "### DONE-B" | tee -a $FRES
