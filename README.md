# mm2-geo — reproducibility (companion repository)

This repository regenerates every table and figure of the mm2-geo paper. The **tool** itself
(minimap2 v2.30 + the geometric patch, with the gated diagnostic flags) lives in a separate
repository; this one holds the scripts, small data artifacts and raw outputs.

- Tool repository: https://github.com/genpat-it/mm2-geo (tag `submission-2026`)
- Manuscript: "Geometry-guided band-doubling: a lightweight heuristic for faster minimap2 gap
  filling with near-identical small-variant calling accuracy".
- Archival: a permanent Zenodo snapshot (this repo + the tool repo) with a DOI will be deposited
  **upon acceptance** (no pre-acceptance DOI, to avoid a permanent orphan record).

## Layout
```
scripts/                one script per table/figure (see map below)
data/                   small artifacts (regenerable / fetched by accession otherwise):
  ecoli_strain_assembly.fasta   Flye assembly of SRR9900640 (spike-in base)
  ecoli_ref_mut.fa              assembly + implanted variants (spike-in reference)
  ecoli_spikein_truth.vcf.gz    planted truth (regenerable via scripts/15, seed 20260716)
  genomewide_readnames.txt      read-name manifest of the 138,688-read genome-wide subset
environment.yml         conda env pinning tool versions (DeepVariant via docker; see file)
ACCESSIONS.txt          all data accessions / sources (raw reads fetched by accession, not bundled)
CHECKSUMS.md5           md5 of the data/ artifacts
geo.patch               copy of the clean geometric patch (canonical copy in the tool repo)
```
The DP-cell measurement (Table 3) and the geometry-vs-doubling ablation (Table 4) require the tool
repo's **gated instrumentation** (`MM2_GEO_STATS`, `MM2_GEO_FIXED`, `MM2_GEO_NODOUBLE`); build the
tool from the mm2-geo repository, which includes them (zero cost when unset). `geo.patch` is the clean
upstreamable change only and does not contain those flags.

## Setup (do this once)
The pipeline uses standard bioinformatics tools; you set up the environment on your machine.

1. **Create the environment** (all analysis tools):
   ```bash
   conda env create -f environment.yml && conda activate mm2-geo-repro
   ```
   DeepVariant is not conda-installable — the three-caller cross-check (Table S5) uses
   `google/deepvariant:1.6.1` from its own image, as in the paper. Everything else is in this env.
2. **Build the tool** (one binary is BOTH stock minimap2 and mm2-geo):
   ```bash
   git clone https://github.com/genpat-it/mm2-geo && cd mm2-geo
   git checkout submission-2026 && make && cp minimap2 mm2geo
   ```
3. **Download the data** by accession into a working directory — see `ACCESSIONS.txt` (GIAB HG002
   HiFi/ONT/CLR, *E. coli* `SRR9900640`, references, truth sets). Raw reads are large and not bundled.
4. **Configure** `scripts/config.sh`: set `WORKDIR` (where you put the data), `MM2GEO` (the binary
   from step 2) and the Clair3/RTG conda env names — or override any of them from the environment.
5. **Verify** (~30 s, synthetic, no downloads):
   ```bash
   bash scripts/smoke_test.sh      # confirms build + samtools + config are wired correctly
   ```
6. **Run** any numbered script, e.g. `bash scripts/01_benchmark_grid.sh` (see the map below).

## Configuration (edit once)
All scripts `source scripts/config.sh`, which defines every path/tool in one place: the working
directory holding the fetched data, the built `mm2geo` binary, the Clair3/RTG conda environment
names, the Clair3 model directories and the thread count. **Edit `scripts/config.sh` (or override any
value from the environment, e.g. `WORKDIR=/data bash scripts/01_benchmark_grid.sh`) to match your
clone before running anything.** Nothing in the scripts changes the analysis — `config.sh` only says
where inputs, tools and the binary live. Tools (`samtools`, `bcftools`, `rtg`, ...) resolve on `PATH`
once the environment from `environment.yml` is active.

## Build and gating
```bash
git clone https://github.com/genpat-it/mm2-geo && cd mm2-geo
git checkout submission-2026   # exact version reported in the paper
make
cp minimap2 mm2geo          # the Makefile emits ./minimap2; we run it as ./mm2geo
```
Environment variables (all default off; unset ⇒ stock minimap2 behaviour):

| variable | effect |
|---|---|
| `MM2_GEO=1` | enable geometry-guided band-doubling |
| `MM2_GEO_MARGIN=m` | band margin (default 20) |
| `MM2_GEO_STATS=1` | print DP-cell count + band histogram to stderr (collect at `-t1` for band stats) |
| `MM2_GEO_FIXED=1` | ablation: fixed band `b=m`, ignore geometry |
| `MM2_GEO_NODOUBLE=1` | ablation: single pass, no retry |

With `MM2_GEO=0` the executable's alignment records are byte-identical to stock minimap2 2.30
(headers aside) — verified by `scripts/00_verify_output_identity.sh`.

## Software versions
minimap2 2.30-r1287 · mm2-fast 2.24-r1122 · samtools 1.21 · bcftools 1.23 · Clair3 v2.0.2
(models `hifi`, `r1041_e82_400bps_sup_v500`, and `..._v430_bacteria_finetuned` for *E. coli*) ·
DeepVariant 1.6.1 (docker, PACBIO) · RTG Tools 3.13 · Sniffles2 2.7.2 · truvari 5.4.0 ·
Flye 2.9.6-b1802 · Badread 0.4.1. Compiler: GCC 14.3.1, `-O2 -Wall`.

## Script → table/figure map
Scripts are numbered in run order; the table/figure each produces is listed here (filenames do
not embed table numbers, so they stay valid if the layout changes).

| Script | Produces |
|---|---|
| `00_verify_output_identity.sh` | `MM2_GEO=0` == minimap2 check (Methods) |
| `01_benchmark_grid.sh` | Table 2 (time / memory / F1) + Figure 2 (F1 bars) |
| `02_dp_cell_count.sh` | Table 3 (logical DP cells + band stats; uses `MM2_GEO_STATS`) |
| `03_ablation_cigar_identity.sh` | Table 4, CIGAR-identity columns |
| `04_ablation_variant_f1.sh` | Table 4, F1 columns (fixed-start / geometry-only / mm2-geo) |
| `05_alignment_concordance.sh` | Table 6 (placement / MAPQ / CIGAR / NM concordance) |
| `06_primary_presence.sh` | Supplementary Table S3 (primary-presence counts) |
| `07_sv_calling.sh` | Table 7 (Sniffles2 + truvari vs GIAB Tier1) |
| `08_composition_mm2fast_f1.sh` | Table 8, F1 (mm2-fast / mm2-fast-geo); timings from `01` |
| `09_margin_sweep.sh` | Table 10 (margin sweep m ∈ {5,10,20,40,80}) |
| `10_thread_scaling.sh` | Table 11 + Figure 3 (1–64 threads) |
| `11_genomewide_speed.sh` | §3.12 genome-wide speed / memory |
| `12_genomewide_concordance.sh` | §3.12 genome-wide placement / MAPQ / CIGAR / NM |
| `13_stress_dense_indels.sh` | §3.10 stress test, dense small indels (Badread) |
| `14_stress_large_sv.sh` | §3.10 stress test, large structural variants (Badread) |
| `15_ecoli_plant_variants.py` | helper: implant the spike-in truth (seed 20260716) |
| `16_ecoli_spikein_f1.sh` | Table 2 *E. coli* row + §3.1 bacterial spike-in (uses `15`) |
| `17_indel_length_stratified_f1.sh` | Supplementary Table S10 (F1 by indel length; calls `17_stratify_indel.py`) |

Note: the three-caller cross-check (Table 5, HiFi chr22) uses clair3, bcftools and DeepVariant on
the same BAM (commands in the manuscript Supplementary); the fixed-`-r` comparison (Table 9) reuses
`01_benchmark_grid.sh` with minimap2 `-r 20`.

## Random seeds
- *E. coli* spike-in: seed `20260716` (hard-coded in `ecoli_plant_variants.py`); regenerates
  `ref_mut.fa` + `truth.vcf` (2,605 SNPs + 724 indels) deterministically. Sanity: applying the
  truth VCF to `ref_mut.fa` with `bcftools consensus` reproduces the original assembly exactly.
- Genome-wide subset: `reads_200k.fq` is a whole-genome random sample of the GIAB HG002 PacBio CCS
  data (138,688 reads, 1.78 Gb, ~0.57× of GRCh38), sampled without reference to mapping location.
  The exact sampling command and seed are recorded alongside the read-name manifest in the archive.

## Data provenance (fetched by accession; not bundled)
- Human GIAB HG002: PacBio HiFi (CCS 15/20 kb, GRCh38), ONT R10.4.1 (2023.05 super-accuracy),
  PacBio CLR (MtSinai). Small-variant truth: GIAB HG002 GRCh38 benchmark VCF + high-confidence BED.
  SV truth: GIAB HG002 SV Tier1 v0.6 (GRCh37).
- *E. coli*: real ONT isolate ENA `SRR9900640`; reference K-12 MG1655 `NC_000913.3`.

## Reproducing a result (pattern)
Every script writes a `results_*.txt`. Example (benchmark grid):
```bash
cd scripts && bash 01_benchmark_grid.sh            # -> results with time, RSS, F1 per row
```
The *E. coli* spike-in (real reads + known truth):
```bash
flye --nano-raw SRR9900640.fastq --out-dir asm --threads 32
python 15_ecoli_plant_variants.py                  # -> ref_mut.fa + truth.vcf (seed 20260716)
bash 16_ecoli_spikein_f1.sh                         # map (mm2/geo) -> clair3 bacterial -> vcfeval
```
F1 stratified by indel length (Supplementary Table S10):
```bash
cd scripts && bash 17_indel_length_stratified_f1.sh # map (MM2_GEO=0 vs =1) -> clair3 -> vcfeval,
                                                    # then 17_stratify_indel.py bins by |REF-ALT| length
```

## Caveat on the E. coli absolute F1
The spike-in truth is exact, but the strain reference is an unpolished Flye assembly with residual
errors that generate background calls **equally** for minimap2 and mm2-geo. The absolute F1
(SNV ~0.90, indel ~0.76) is therefore not comparable to the GIAB human values; the **paired
minimap2-vs-mm2-geo difference** (ΔSNV +0.0018, ΔINDEL −0.0087) is the intended measure.
