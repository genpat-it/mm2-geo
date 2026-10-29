# mm2-geo — reproducibility material

This branch regenerates every table and figure of the mm2-geo paper ("Geometry-guided band sizing with a
score certificate: faster minimap2 gap filling that preserves its gap scores"). The tool itself
(minimap2 v2.30 plus a patch to `align.c`) is on the `main` branch of https://github.com/genpat-it/mm2-geo.
The tagged tree also adds a DP-cell counter in `ksw2_extd2_sse.c` for `MM2_GEO_STATS`; `geo.patch` in the tool
repository is the minimal production patch (`align.c` only, without the diagnostic flags).

| tag | content |
|---|---|
| `revision-2026` (tool) | certified mode (`MM2_GEO=1`, the default) and heuristic mode (`MM2_GEO=1 MM2_GEO_FAST=1`) |
| `submission-2026` (tool) | heuristic mode only, as first released; unchanged in `revision-2026` |
| `reproducibility-revision-2026` (this branch) | scripts for every table of the paper and its Supplementary Material |
| `reproducibility-submission-2026` (this branch) | the material of the first submission |

A permanent Zenodo snapshot of both branches will be deposited upon acceptance.

## Requirements

- Linux x86-64; gcc >= 9 (C99), make, zlib headers, git, curl, GNU time (`/usr/bin/time`), bc
- conda (the environment in `environment.yml` provides samtools, bcftools, Clair3, RTG Tools, Sniffles2,
  truvari, Flye, Badread, pbsim3)
- docker for DeepVariant in script 32 (or `DEEPVARIANT=skip`)
- a CPU with AVX2 for mm2-fast (scripts 01, 08, 10, 28; AVX-512 is used when present)
- disk: ~400 GB for everything, ~150 GB for the core scripts (`run_all.sh --core`)

## Quick start (from an empty machine)

```bash
git clone -b reproducibility https://github.com/genpat-it/mm2-geo mm2-geo-repro && cd mm2-geo-repro
conda env create -f environment.yml && conda activate mm2-geo-repro
export WORKDIR=/path/with/space          # ~150 GB for "core", ~400 GB with the genome-wide data
bash scripts/build_tools.sh              # builds every executable from the public tags into $WORKDIR/bin
bash scripts/fetch_data.sh core          # references, truth sets, stratifications, chr22/chr14 read slices
bash scripts/smoke_test.sh               # ~30 s sanity check of the build and configuration
bash scripts/00_verify_output_identity.sh
# ... then any numbered script; the genome-wide ones need:
bash scripts/fetch_data.sh genomewide    # complete GIAB HG002 HiFi BAM (121 GB) and ONT R10.4.1 CRAM (104 GB)
```

Or everything in one command (build, data, smoke test and every numbered script in dependency order; one log
per step and a summary of exit codes in `$WORKDIR/logs/status.tsv`; completed steps are skipped on a re-run):

```bash
bash scripts/run_all.sh            # all tables (~2 days, ~400 GB); --core skips the genome-wide scripts 24, 27, 34 and 35
```

Everything the scripts need (paths, executables, conda environments, Clair3 models, threads) is defined
once in `scripts/config.sh`; any value can be overridden from the environment. Results are written to
`$WORKDIR/hifi_bench/` (scripts 01–17, 32 and 33) and `$WORKDIR/results/` (scripts 00, 18–31 and 34; e.g.
`results/cert_oracle/` for 30 and 31, whose table also goes to stdout, `results/ont_r10_genomewide/` for 27 and
`results/ont_r9_genomewide/` for 34). Timings assume an otherwise idle machine and local storage for the inputs; all
timing scripts use one executable per comparison, change only environment variables, discard the SAM
output and report the median of three interleaved runs.

## Executables (`build_tools.sh`)

| config variable | built from | used for |
|---|---|---|
| `MM2GEO` | tag `revision-2026` | all heuristic- and certified-mode results |
| `MM2GEO_SUB` | tag `submission-2026` | heuristic-mode speed checks with the code as first released |
| `MM2GEO_INSTR` | `revision-2026` + `instrumentation.patch` | per-gap log (`MM2_GEO_GAPLOG`) and DP-stage CPU time (`MM2_GEO_TIMEDP`); output unchanged |
| `MM2GEO_ORACLE` | `revision-2026` + `cert_oracle.patch` | per-gap check of the certified mode against the full band (`MM2_GEO_CERT_ORACLE`); output unchanged |
| `MM2GEO_O3` | `revision-2026` with bioconda's minimap2 compiler flags | build comparison |
| `MM2FAST`, `MM2FASTGEO` | mm2-fast commit `14fe36c`, without / with `geo_mm2fast.patch` (AVX-512 build when the CPU supports it, else AVX2) | composition with mm2-fast (optional) |

Environment variables of the tool (all off by default; unset means stock minimap2):

| variable | effect |
|---|---|
| `MM2_GEO=1` | certified mode (default): geometry-sized band accepted by the score certificate, exact one-gap shortcut, score-only first pass; output identical to minimap2 |
| `MM2_GEO_FAST=1` | with `MM2_GEO=1`: heuristic mode (geometry-sized band, doubled on boundary contact), as selected by `MM2_GEO=1` alone in `submission-2026` |
| `MM2_GEO_CERT=1` | accepted for compatibility; no effect |
| `MM2_GEO_MARGIN=m` | band margin (default 20) |
| `MM2_GEO_ONEGAP=1` | with `MM2_GEO_FAST=1`: adds the exact one-gap shortcut alone |
| `MM2_GEO_CERT_SCOREPROBE=n` | certified mode: minimum gap length for the score-only first pass (default 256; 0 disables) |
| `MM2_GEO_CERT_ORACLE=1` / `=2` | oracle build only (`cert_oracle.patch`): per-gap check of the certified / heuristic mode against the full band |
| `MM2_GEO_STATS=1`, `MM2_GEO_FIXED=1`, `MM2_GEO_NODOUBLE=1` | DP-cell statistics and ablations |
| `MM2_GEO_GAPLOG=file`, `MM2_GEO_TIMEDP=1` | instrumented build only: per-gap log, DP-stage CPU time |

## Script → table map

Table numbers refer to the revised manuscript (main text: Tables 1–3, Figure 1; Supplementary: S1–S23).

| script | produces |
|---|---|
| `00_verify_output_identity.sh` | output identity: unset = minimap2 v2.30, certified = minimap2, heuristic unchanged across tags |
| `01_benchmark_grid.sh` | Table 1 (heuristic mode: time, memory, small-variant F1); Table S13 timing rows (mm2-fast, mm2-fast+geo) |
| `02_dp_cell_count.sh` | Table S2 (logical DP cells) |
| `03_ablation_cigar_identity.sh`, `04_ablation_variant_f1.sh` | Table 3 (geometry vs doubling ablation) |
| `05_alignment_concordance.sh`, `06_primary_presence.sh` | Tables S6, S4 |
| `07_sv_calling.sh` | Table S7 (SVs vs GIAB Tier1, heuristic mode) |
| `08_composition_mm2fast_f1.sh` | Table S13 (F1 with mm2-fast) |
| `09_margin_sweep.sh` | Table S9 |
| `10_thread_scaling.sh` | Table S12 |
| `11_genomewide_speed.sh`, `12_genomewide_concordance.sh` | Table S14 (0.57× genome-wide subset) |
| `13_stress_dense_indels.sh`, `14_stress_large_sv.sh` | Supplementary section "Stress tests" |
| `15_ecoli_plant_variants.py`, `16_ecoli_spikein_f1.sh` | *E. coli* spike-in (Table S3 and Supplementary text) |
| `17_indel_length_stratified_f1.sh` | Table S10 |
| `18_gap_rounds_dE.sh` | Table S15 (band-doubling rounds, joint (d,E) statistics) |
| `19_discordant_reads.sh` | Table S16 (discordant reads: severity, GIAB stratifications); uses the output of 18 |
| `20_simulated_truth.sh` | discordant reads against simulated truth (pbsim3, seed 2024) |
| `21_callset_diff.sh` | Table S17 (SV and small-variant callsets on the same GRCh38 alignments) |
| `22_sv_truth_grch37.sh` | Table S17, bottom block (SVs vs GIAB Tier1 for minimap2 / heuristic / certified) |
| `23_modes_identity_timing.sh` | Table 2 end-to-end columns, Table S22 (both modes, certified identity) |
| `24_genomewide_perchr.sh` | Table S18 (full-coverage HG002 HiFi per chromosome, both modes; chr1-only setting) |
| `25_speed_threads_builds.sh` | Table S20 (threads, output handling, index scope, bioconda and -O3 builds) |
| `26_dp_stage.sh` | Table 2 DP columns, Table S21 (DP-stage decomposition, Amdahl prediction) |
| `27_ont_r10_genomewide.sh` | ONT R10.4.1 genome-wide subset (Tables S19, S21, S22) |
| `28_mm2fast_modes.sh` | Table S22, mm2-fast rows |
| `29_input_counts.sh` | Table S19 (input reads, bases, alignment rates) |
| `30_cert_oracle.sh` | Table S23 (per-gap oracle of the certified mode against the full band; heuristic-mode positive control) |
| `31_bound_exhaustive.py` | Supplementary proof section, numerical check (exhaustive DP over all paths of small gaps, three scorings) |
| `32_three_callers.sh` | Table S5 (Clair3, bcftools and DeepVariant on minimap2 and heuristic-mode BAMs; DeepVariant via docker) |
| `33_fixed_r_band.sh` | Table S8 (minimap2 `-r 20,20` vs default vs heuristic mode: time and F1) |
| `34_ont_r9_genomewide.sh` | ONT R9.4.1 genome-wide subset (Tables S19, S21, S22; first 20 GB of the GIAB ultra-long Guppy 3.2.4 FASTQ) |
| `35_hifi53x_alignment_rate.sh` | Table S19, full-coverage HiFi row: alignment rate of stock minimap2 on all 11.7 million reads |

Tables S1 (environment variables) and S11 (aligners compared) are descriptive; Figure 1 is a schematic.

## Raw outputs
The summary files of the runs reported in the paper are in `outputs/` (see `outputs/README.md` for the
file → table map); a re-run writes the same files into `$WORKDIR`.

## Data
`scripts/fetch_data.sh` downloads and prepares every input from public sources (listed in
`ACCESSIONS.txt`): GRCh38 (GenBank `GCA_000001405.15`), the GRCh38 no-alt analysis set (to decode the ONT
CRAM locally), GRCh37 chr22, GIAB v4.2.1 small-variant truth, GIAB SV Tier1 v0.6, GIAB v3.3
stratifications, chr22/chr14 slices of the GIAB HG002 HiFi, ONT R10.4.1 and CLR alignments (streamed),
the genome-wide HiFi set (one complete SMRT Cell of the GIAB HG002 CCS 15 kb data, 138,688 reads, ~0.57×, md5-checked), the *E. coli* ONT isolate
(ENA `SRR9900640`, md5-checked) and K-12 MG1655 reference (`NC_000913.3`), pbsim3 reads simulated from
chr22 (seed 2024), and, with `genomewide`, the complete HiFi BAM and ONT CRAM, and the ONT R9.4.1 genome-wide subset (streamed from GIAB, checked by read and base counts). Small artifacts are in
`data/` (*E. coli* strain assembly, spike-in reference and truth; md5 in
`CHECKSUMS.md5`). `build_tools.sh` builds mm2-geo from the tagged source tree itself (not by applying
`geo.patch`), and takes `geo_mm2fast.patch` from tag `revision-2026` to build mm2-fast+geo, so there is a single
copy of each patch (in the tool repository); `instrumentation.patch` and `cert_oracle.patch` (analysis builds)
are in this branch.

## Software
minimap2 2.30-r1287 (`v2.30` of lh3/minimap2) · mm2-fast 2.24-r1122 (commit `14fe36c`) · samtools 1.21 ·
bcftools 1.23 · Clair3 v2.0.2 (models `hifi`, `r1041_e82_400bps_sup_v500`) · DeepVariant 1.6.1 · RTG Tools
3.13 · Sniffles2 2.7.2 · truvari 5.4.0 · Flye 2.9.6 · Badread 0.4.1 · pbsim3 3.0 · GCC 14.3.1 (`-O2 -Wall`).

## Seeds
*E. coli* spike-in `20260716`; simulated chr22 reads `2024`; stress tests (scripts 13 and 14: implanted
variants and Badread reads) `77`; ONT R10.4.1 genome-wide subset `samtools -s 42.10` (seed 42, 10%).
