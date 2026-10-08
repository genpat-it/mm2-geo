#!/bin/bash
# Variant callsets on the same GRCh38 chr22 alignments (Supplementary Table S17): SVs (Sniffles2) and small
# variants (Clair3, whole chr22, no BED), heuristic and certified modes compared with the minimap2 callset
# (truvari / rtg vcfeval with minimap2 as baseline), plus the high-confidence BED restriction.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; D=results/callset; mkdir -p $D
source "$CONDA_PROFILE"
for s in "hifi map-hifi hifi_chr22.fq hifi $MODEL_HIFI" "ont map-ont ont_r10_chr22.fq ont $MODEL_ONT" "clr map-pb clr_chr22.fq none none"; do set -- $s
  for c in stock heuristic certified; do
    case $c in stock) E="MM2_GEO=0";; heuristic) E="MM2_GEO=1";; certified) E="MM2_GEO=1 MM2_GEO_CERT=1";; esac
    b=$D/${1}_$c.bam
    [ -s $b.bai ] || { env $E "$MM2GEO" -ax $2 -t"$THREADS" chr22_named.fa $3 2>/dev/null | "$SAMTOOLS" sort -@8 -m2G -o $b - 2>/dev/null
      "$SAMTOOLS" reheader -c 'sed "s/SN:chr22/SN:CM000684.2/"' $b > $b.tmp && mv $b.tmp $b; "$SAMTOOLS" index $b; }
    [ -s $D/${1}_$c.sv.vcf.gz ] || "$SNIFFLES" --input $b --vcf $D/${1}_$c.sv.vcf.gz --reference chr22_cm.fa --threads 16 > $D/snf_${1}_$c.log 2>&1
    if [ $4 != none ] && [ ! -s $D/c3_${1}_$c/merge_output.vcf.gz ]; then conda activate "$CONDA_ENV_CLAIR3"
      "$CLAIR3" --bam_fn=$b --ref_fn=chr22_cm.fa --threads="$THREADS" --platform=$4 --model_path=$5 --output=$D/c3_${1}_$c --ctg_name=CM000684.2 > $D/c3_${1}_$c.log 2>&1
      conda deactivate; fi
  done
  echo "== $1: BAM records identical to minimap2 (certified): $(cmp -s <("$SAMTOOLS" view $D/${1}_stock.bam) <("$SAMTOOLS" view $D/${1}_certified.bam) && echo yes || echo no)" | tee -a $D/summary.txt
  for c in heuristic certified; do
    rm -rf $D/sv_${1}_$c; "$TRUVARI" bench -b $D/${1}_stock.sv.vcf.gz -c $D/${1}_$c.sv.vcf.gz -o $D/sv_${1}_$c --passonly > /dev/null 2>&1
    python3 -c "import json;d=json.load(open('$D/sv_${1}_$c/summary.json'));print('   SV $c: only-minimap2 %d only-geo %d shared %d'%(d['FN'],d['FP'],d['TP-base']))" | tee -a $D/summary.txt
    if [ $4 != none ]; then conda activate "$CONDA_ENV_RTG"
      for x in all hc; do bed=""; [ $x = hc ] && bed="-e conf_chr22.bed"; rm -rf $D/sm_${x}_${1}_$c
        "$RTG" vcfeval -b $D/c3_${1}_stock/merge_output.vcf.gz -c $D/c3_${1}_$c/merge_output.vcf.gz -t ref_chr22_sdf $bed -o $D/sm_${x}_${1}_$c > /dev/null 2>&1
        awk -v c=$c -v x=$x '$1=="None"{printf "   small variants %s (%s): only-minimap2 %s only-geo %s shared %s\n", c, x, $5, $4, $3}' $D/sm_${x}_${1}_$c/summary.txt | tee -a $D/summary.txt; done
      conda deactivate; fi
  done
done
