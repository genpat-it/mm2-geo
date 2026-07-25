# config.sh — EDIT THIS FILE, then run any numbered script; they all `source` it.
#
# Every path/tool the reproduction scripts need is defined here once. Values can also be
# overridden from the environment (e.g. `WORKDIR=/data bash 01_benchmark_grid.sh`).
# Nothing below changes the analysis — only where inputs, tools and the binary are found.

# ---------------------------------------------------------------------------
# 1) Working directory: where the reference/reads/truth artifacts live.
#    Fetch the data per ACCESSIONS.txt into this directory (chr22_named.fa,
#    chr22_cm.fa, truth_chr22.vcf.gz, conf_chr22.bed, ref_chr22_sdf, hifi_chr22.fq,
#    ont_r10_chr22.fq, ...). Scripts `cd` here, so those names are relative to it.
WORKDIR="${WORKDIR:-/path/to/mm2geo-workdir}"   # <-- SET THIS (or export WORKDIR=...) before running

# ---------------------------------------------------------------------------
# 2) The mm2-geo binary, built from github.com/genpat-it/mm2-geo (tag submission-2026):
#      git clone ... && cd mm2-geo && git checkout submission-2026 && make && cp minimap2 mm2geo
#    Point MM2GEO at that executable (one binary; MM2_GEO=0 == stock minimap2).
MM2GEO="${MM2GEO:-$WORKDIR/mm2geo}"
# The mm2-geo source tree (only script 00 needs it: it rebuilds with the patch reverted to
# prove MM2_GEO=0 is byte-identical to stock minimap2).
MM2GEO_SRC="${MM2GEO_SRC:-$WORKDIR/../mm2geo_src}"
# Optional: mm2-fast and mm2-fast+geo binaries (scripts 01/10 only; build from
# bwa-mem2/mm2-fast v2.24 with geo.patch applied). Leave as-is if not reproducing those rows.
MM2FAST="${MM2FAST:-$WORKDIR/../mm2-fast/mm2fast}"
MM2FASTGEO="${MM2FASTGEO:-$WORKDIR/../mm2-fast/mm2fastgeo}"

# ---------------------------------------------------------------------------
# 3) Tools. With the conda env from environment.yml active they are on PATH, so the
#    bare names work; otherwise set absolute paths here.
SAMTOOLS="${SAMTOOLS:-samtools}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
RTG="${RTG:-rtg}"
BGZIP="${BGZIP:-bgzip}"
TABIX="${TABIX:-tabix}"
SNIFFLES="${SNIFFLES:-sniffles}"   # structural-variant caller (script 07)
TRUVARI="${TRUVARI:-truvari}"      # SV benchmarking (script 07)

# ---------------------------------------------------------------------------
# 4) Conda: profile to source + environment names holding Clair3 and RTG Tools.
#    (environment.yml ships a single env; set both names to it, or create separate envs.)
CONDA_PROFILE="${CONDA_PROFILE:-$HOME/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV_CLAIR3="${CONDA_ENV_CLAIR3:-c3}"
CONDA_ENV_RTG="${CONDA_ENV_RTG:-rtg}"
CLAIR3="${CLAIR3:-run_clair3.sh}"   # resolves on PATH once CONDA_ENV_CLAIR3 is active

# ---------------------------------------------------------------------------
# 5) Clair3 model directories (ship with Clair3; adjust to your install).
CLAIR3_MODELS="${CLAIR3_MODELS:-$HOME/miniconda3/envs/c3/bin/models}"
MODEL_HIFI="${MODEL_HIFI:-$CLAIR3_MODELS/hifi}"
MODEL_ONT="${MODEL_ONT:-$CLAIR3_MODELS/r1041_e82_400bps_sup_v500}"
MODEL_BACT="${MODEL_BACT:-$CLAIR3_MODELS/r1041_e82_400bps_sup_v430_bacteria_finetuned}"
# E. coli K-12 MG1655 reference (NC_000913.3), used only for the E. coli concordance row.
ECOLI_REF="${ECOLI_REF:-$WORKDIR/ecoli_MG1655.fna}"
# Full GRCh38 + genome-wide HG002 HiFi subset (scripts 11/12).
GRCH38_REF="${GRCH38_REF:-$WORKDIR/GRCh38.fa}"
GW_READS="${GW_READS:-$WORKDIR/reads_200k.fq}"

# ---------------------------------------------------------------------------
# 6) Misc.
THREADS="${THREADS:-32}"
export TMPDIR="${TMPDIR:-$WORKDIR/tmp}"
mkdir -p "$TMPDIR"
