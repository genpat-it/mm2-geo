# config.sh — EDIT THIS FILE, then run any numbered script; they all `source` it.
#
# Every path/tool the reproduction scripts need is defined here once. Values can also be
# overridden from the environment (e.g. `WORKDIR=/data bash 01_benchmark_grid.sh`).
# Nothing below changes the analysis — only where inputs, tools and the binary are found.

# Absolute paths of this repository (resolved here, before any script changes directory).
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(dirname "$SCRIPTS")"

# ---------------------------------------------------------------------------
# 1) Working directory: where the reference/reads/truth artifacts live.
#    `bash scripts/fetch_data.sh` downloads and prepares every input here (chr22_named.fa,
#    chr22_cm.fa, truth_chr22.vcf.gz, conf_chr22.bed, ref_chr22_sdf, hifi_chr22.fq,
#    ont_r10_chr22.fq, ...; see ACCESSIONS.txt). Scripts `cd` here, so those names are relative to it.
WORKDIR="${WORKDIR:-/path/to/mm2geo-workdir}"   # <-- SET THIS (or export WORKDIR=...) before running

# ---------------------------------------------------------------------------
# 2) Executables. `bash scripts/build_tools.sh` builds all of them into $WORKDIR/bin from the tool
#    repository (github.com/genpat-it/mm2-geo); one executable is both stock minimap2 (MM2_GEO unset)
#    and mm2-geo.
#      MM2GEO        tag revision-2026: heuristic mode (MM2_GEO=1) and certified mode (MM2_GEO=1 MM2_GEO_CERT=1)
#      MM2GEO_SUB    tag submission-2026: the originally released heuristic mode (unchanged
#                    in revision-2026); the heuristic-mode speed scripts can be run with either
#      MM2GEO_INSTR  revision-2026 + instrumentation.patch: adds MM2_GEO_GAPLOG (per-gap log) and
#                    MM2_GEO_TIMEDP (DP-stage CPU time); analysis only, same output
#      MM2GEO_O3     revision-2026 rebuilt with bioconda's minimap2 compiler flags (-O3 ...)
BIN="${BIN:-$WORKDIR/bin}"
MM2GEO="${MM2GEO:-$BIN/mm2geo}"
MM2GEO_SUB="${MM2GEO_SUB:-$BIN/mm2geo_submission}"
MM2GEO_INSTR="${MM2GEO_INSTR:-$BIN/mm2geo_instr}"
MM2GEO_O3="${MM2GEO_O3:-$BIN/mm2geo_o3}"
# stock minimap2 2.30 from bioconda (build-comparison script only): `conda install -c bioconda minimap2=2.30`
BIOCONDA_MM2="${BIOCONDA_MM2:-minimap2}"
MM2GEO_SRC="${MM2GEO_SRC:-$BIN/src/mm2-geo}"
# Optional: mm2-fast and mm2-fast+geo binaries (scripts 01, 08, 10, 28), built by build_tools.sh from
# bwa-mem2/mm2-fast commit 14fe36c without / with geo_mm2fast.patch.
MM2FAST="${MM2FAST:-$BIN/mm2fast}"
MM2FASTGEO="${MM2FASTGEO:-$BIN/mm2fastgeo}"   # mm2-fast + geo_mm2fast.patch (heuristic and certified modes)

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
# 4) Conda: profile to source + environments holding Clair3 and RTG Tools. By default both are the
#    environment that is active when a script starts (the one created from environment.yml).
_CONDA_BASE="${CONDA_EXE:+$(dirname "$(dirname "$CONDA_EXE")")}"; _CONDA_BASE="${_CONDA_BASE:-$HOME/miniconda3}"
CONDA_PROFILE="${CONDA_PROFILE:-$_CONDA_BASE/etc/profile.d/conda.sh}"
CONDA_ENV_CLAIR3="${CONDA_ENV_CLAIR3:-${CONDA_PREFIX:-mm2-geo-repro}}"
CONDA_ENV_RTG="${CONDA_ENV_RTG:-${CONDA_PREFIX:-mm2-geo-repro}}"
CLAIR3="${CLAIR3:-run_clair3.sh}"   # resolves on PATH once CONDA_ENV_CLAIR3 is active
# conda_on ENV / conda_off: switch to ENV for one tool and back. When ENV is already the active environment
# (the default, a single environment from environment.yml) both are no-ops, so that a plain `conda deactivate`
# can never drop the scripts into the base environment.
conda_on(){ _CONDA_SWITCHED=0
  case "$1" in "${CONDA_PREFIX:-}"|"${CONDA_DEFAULT_ENV:-}") return 0;; esac
  type conda >/dev/null 2>&1 || source "$CONDA_PROFILE"
  conda activate "$1" && _CONDA_SWITCHED=1; }
conda_off(){ [ "${_CONDA_SWITCHED:-0}" = 1 ] && conda deactivate; _CONDA_SWITCHED=0; }

# ---------------------------------------------------------------------------
# 5) Clair3 model directories (shipped with the bioconda Clair3 package in <env>/bin/models).
CLAIR3_MODELS="${CLAIR3_MODELS:-${CONDA_PREFIX:-$_CONDA_BASE/envs/mm2-geo-repro}/bin/models}"
MODEL_HIFI="${MODEL_HIFI:-$CLAIR3_MODELS/hifi}"
MODEL_ONT="${MODEL_ONT:-$CLAIR3_MODELS/r1041_e82_400bps_sup_v500}"
MODEL_BACT="${MODEL_BACT:-$CLAIR3_MODELS/r1041_e82_400bps_sup_v430_bacteria_finetuned}"
# E. coli K-12 MG1655 reference (NC_000913.3), used only for the E. coli concordance row.
ECOLI_REF="${ECOLI_REF:-$WORKDIR/ecoli_MG1655.fna}"
# Full GRCh38 + genome-wide HG002 HiFi subset (scripts 11/12).
GRCH38_REF="${GRCH38_REF:-$WORKDIR/GRCh38.fa}"
GW_READS="${GW_READS:-$WORKDIR/hifi_genomewide_subset.fq}"   # GIAB HG002 CCS 15 kb, m54238_180901_011437.Q20.fastq
# Complete GIAB HG002 data sets for the genome-wide benchmarks (fetched by `fetch_data.sh genomewide`).
HG002_HIFI_BAM="${HG002_HIFI_BAM:-$WORKDIR/HG002.SequelII.merged_15kb_20kb.GRCh38.duplomap.bam}"
HG002_ONT_CRAM="${HG002_ONT_CRAM:-$WORKDIR/PAO83395.pass.cram}"
# pbsim3 reads simulated from GRCh38 chr22 with their true alignments (MAF), for the simulated-truth check.
SIM_READS="${SIM_READS:-$WORKDIR/sim/chr22sim.fq}"
SIM_MAF="${SIM_MAF:-$WORKDIR/sim/chr22sim.maf}"
PBSIM="${PBSIM:-pbsim}"   # pbsim3 v3.0 (bioconda: pbsim3)
PBSIM_MODEL="${PBSIM_MODEL:-$CONDA_PREFIX/data/ERRHMM-SEQUEL.model}"

# ---------------------------------------------------------------------------
# 6) Misc.
THREADS="${THREADS:-32}"
export TMPDIR="${TMPDIR:-$WORKDIR/tmp}"
# scripts 01-17 write their tables to $WORKDIR/hifi_bench, scripts 00 and 18-29 to $WORKDIR/results
mkdir -p "$TMPDIR" "$WORKDIR/hifi_bench" "$WORKDIR/results"
