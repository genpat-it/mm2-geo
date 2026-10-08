#!/bin/bash
# Fetch and prepare every input used by the numbered scripts into $WORKDIR (see config.sh).
# Usage: bash scripts/fetch_data.sh [core|genomewide|all]
#   core       (default) references, truth sets, stratifications, the chr22/chr14 HiFi/ONT/CLR read
#              slices, the 0.57x genome-wide HiFi subset, the E. coli data and the simulated chr22 reads
#              (scripts 00-22 and 29; the slices stream from the GIAB/ONT servers, so only the selected
#              chromosomes are transferred)
#   genomewide additionally the complete GIAB HG002 HiFi BAM (121 GB) and ONT R10.4.1 CRAM (104 GB) used
#              by the genome-wide scripts (24, 27); needs ~250 GB of free space
#   all        both
# Every step is skipped when its output already exists, so the script can be re-run after an interruption.
set -euo pipefail
source "$(dirname "$0")/config.sh"
MODE="${1:-core}"
mkdir -p "$WORKDIR" "$TMPDIR"; cd "$WORKDIR"
S="$SAMTOOLS"
log(){ echo "[fetch $(date +%T)] $*"; }
get(){ [ -s "$2" ] || { log "download $2"; curl -fL --retry 5 -C - -o "$2.part" "$1" && mv "$2.part" "$2"; }; }

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
HIFI_BAM_URL=$GIAB/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.SequelII.merged_15kb_20kb.GRCh38.duplomap.bam
CLR_BAM_URL=$GIAB/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_minimap2_bam/HG002_PacBio_GRCh38.bam
ONT_CRAM_URL=https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/hg002/sup/PAO83395.pass.cram
TRUTH_URL=$GIAB/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
SV_URL=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6
STRAT_URL=$GIAB/release/genome-stratifications/v3.3/GRCh38@all
GRCH38_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz
GRCH37_CHR22_URL=http://ftp.ensembl.org/pub/grch37/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz
NOALT_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# ---------- references ----------
if [ ! -s GRCh38.fa ]; then get "$GRCH38_URL" GRCh38.fna.gz; log "GRCh38.fa"; gunzip -c GRCh38.fna.gz > GRCh38.fa; fi
[ -s GRCh38.fa.fai ] || $S faidx GRCh38.fa
[ -s chr22_cm.fa ]    || $S faidx GRCh38.fa CM000684.2 > chr22_cm.fa
[ -s chr22_named.fa ] || sed 's/^>CM000684.2.*/>chr22/' chr22_cm.fa > chr22_named.fa
[ -s ref_chr14.fa ]   || $S faidx GRCh38.fa CM000676.2 > ref_chr14.fa
# the ONT CRAM is encoded against the GRCh38 no-alt analysis set; decode it locally with this reference
if [ ! -s GRCh38_no_alt_analysis_set.fna ]; then get "$NOALT_URL" noalt.fna.gz; gunzip -c noalt.fna.gz > GRCh38_no_alt_analysis_set.fna; fi
[ -s GRCh38_no_alt_analysis_set.fna.fai ] || $S faidx GRCh38_no_alt_analysis_set.fna
for f in chr22_cm.fa chr22_named.fa ref_chr14.fa; do [ -s $f.fai ] || $S faidx $f; done
[ -d ref_chr22_sdf ] || "$RTG" format -o ref_chr22_sdf chr22_cm.fa
[ -d ref_chr14_sdf ] || "$RTG" format -o ref_chr14_sdf ref_chr14.fa

# ---------- small-variant truth (GIAB v4.2.1), renamed to the GenBank contig names ----------
get "$TRUTH_URL/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
get "$TRUTH_URL/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
get "$TRUTH_URL/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
for c in "22 CM000684.2" "14 CM000676.2"; do set -- $c
  [ -s truth_chr$1.vcf.gz ] || { echo -e "chr$1\t$2" > rn_chr$1.txt
    "$BCFTOOLS" view -r chr$1 HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | "$BCFTOOLS" annotate --rename-chrs rn_chr$1.txt -Oz -o truth_chr$1.vcf.gz
    "$TABIX" -f -p vcf truth_chr$1.vcf.gz; }
  [ -s conf_chr$1.bed ] || awk -v c=chr$1 -v n=$2 'BEGIN{OFS="\t"} $1==c{$1=n; print}' HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > conf_chr$1.bed
done

# ---------- GIAB v3.3 stratifications (discordant-read analysis) ----------
mkdir -p strat
get "$STRAT_URL/SegmentalDuplications/GRCh38_segdups.bed.gz" strat/GRCh38_segdups.bed.gz
get "$STRAT_URL/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz" strat/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz
get "$STRAT_URL/Mappability/GRCh38_lowmappabilityall.bed.gz" strat/GRCh38_lowmappabilityall.bed.gz

# ---------- SV truth (GIAB Tier1 v0.6, GRCh37) + GRCh37 chr22 ----------
mkdir -p sv
get "$SV_URL/HG002_SVs_Tier1_v0.6.vcf.gz" sv/HG002_SVs_Tier1_v0.6.vcf.gz
get "$SV_URL/HG002_SVs_Tier1_v0.6.vcf.gz.tbi" sv/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
get "$SV_URL/HG002_SVs_Tier1_v0.6.bed" sv/HG002_SVs_Tier1_v0.6.bed
[ -s sv/grch37_chr22.fa ] || { get "$GRCH37_CHR22_URL" sv/grch37_chr22.fa.gz; gunzip -c sv/grch37_chr22.fa.gz > sv/grch37_chr22.fa; $S faidx sv/grch37_chr22.fa; }

# ---------- chromosome read slices (streamed from the remote alignments) ----------
slice(){ # $1=url $2=region $3=out.fq
  [ -s "$3" ] || { log "slice $2 -> $3"; $S view -@8 -b -T GRCh38_no_alt_analysis_set.fna -F 0x900 "$1" "$2" | $S fastq -@8 - > "$3.part" && mv "$3.part" "$3"; }; }
slice "$HIFI_BAM_URL" chr22 hifi_chr22.fq
slice "$HIFI_BAM_URL" chr14 hifi_chr14.fq
slice "$ONT_CRAM_URL" chr22 ont_r10_chr22.fq
slice "$ONT_CRAM_URL" chr14 ont_r10_chr14.fq
slice "$CLR_BAM_URL"  chr22 clr_chr22.fq
slice "$CLR_BAM_URL"  chr14 clr_chr14.fq

# ---------- genome-wide 0.57x HiFi subset (read-name manifest, used by scripts 11/12/26) ----------
if [ ! -s "$GW_READS" ]; then
  log "genome-wide subset from the read-name manifest (streams the full HiFi BAM once)"
  MAN="$(dirname "$0")/../data/genomewide_readnames.txt"
  $S view -@8 -b -F 0x900 -N "$MAN" "$HIFI_BAM_URL" | $S fastq -@8 - > "$GW_READS.part" && mv "$GW_READS.part" "$GW_READS"
fi

# ---------- E. coli: real ONT isolate (ENA SRR9900640), K-12 MG1655 reference, spike-in artifacts ----------
# (scripts 06, 15, 16). The spike-in reference/truth and the strain assembly they derive from are in data/.
DATA="$(cd "$(dirname "$0")/../data" && pwd)"
mkdir -p ecoli_real2/flye ecoli_vc
if [ ! -s ecoli_real2/SRR9900640.fastq ]; then
  get https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR990/000/SRR9900640/SRR9900640_1.fastq.gz ecoli_real2/SRR9900640_1.fastq.gz
  echo "733b61ae7289db93cc46e14dc473122d  ecoli_real2/SRR9900640_1.fastq.gz" | md5sum -c --quiet
  gunzip -c ecoli_real2/SRR9900640_1.fastq.gz > ecoli_real2/SRR9900640.fastq.part && mv ecoli_real2/SRR9900640.fastq.part ecoli_real2/SRR9900640.fastq
fi
get "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text" "$ECOLI_REF"
[ -s ecoli_real2/flye/assembly.fasta ] || cp "$DATA/ecoli_strain_assembly.fasta" ecoli_real2/flye/assembly.fasta
[ -s ecoli_vc/ref_mut.fa ]            || cp "$DATA/ecoli_ref_mut.fa" ecoli_vc/ref_mut.fa
[ -s ecoli_vc/truth.vcf.gz ]          || { cp "$DATA/ecoli_spikein_truth.vcf.gz" ecoli_vc/truth.vcf.gz; "$TABIX" -f -p vcf ecoli_vc/truth.vcf.gz; }

# ---------- complete data sets for the genome-wide benchmarks ----------
if [ "$MODE" = genomewide ] || [ "$MODE" = all ]; then
  get "$HIFI_BAM_URL" "$HG002_HIFI_BAM"; get "$HIFI_BAM_URL.bai" "$HG002_HIFI_BAM.bai"
  get "$ONT_CRAM_URL" "$HG002_ONT_CRAM"; get "$ONT_CRAM_URL.crai" "$HG002_ONT_CRAM.crai"
fi
# ---------- simulated reads with known true alignments (pbsim3 v3.0, chr22, seed 2024) ----------
if [ ! -s "$SIM_READS" ] && command -v "$PBSIM" >/dev/null 2>&1; then
  mkdir -p "$(dirname "$SIM_READS")"; ( cd "$(dirname "$SIM_READS")"
    "$PBSIM" --strategy wgs --method errhmm --errhmm "$PBSIM_MODEL" --genome "$WORKDIR/chr22_named.fa" \
      --depth 5 --length-mean 15000 --seed 2024 --prefix chr22sim > pbsim.log 2>&1
    gunzip -c chr22sim_0001.fq.gz > "$SIM_READS"; gunzip -c chr22sim_0001.maf.gz > "$SIM_MAF" )
fi
log "done ($MODE)"
