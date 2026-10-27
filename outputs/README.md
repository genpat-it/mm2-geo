# Raw outputs

Summary files as produced by the runs reported in the paper. Large intermediates (BAM, SAM/PAF, read files)
are not included; the scripts in `scripts/` regenerate them. In these files `fast` and `geo` denote the
heuristic mode (`MM2_GEO=1 MM2_GEO_FAST=1`; `MM2_GEO=1` in `submission-2026`), `cert` the certified mode (`MM2_GEO=1`, the default in `revision-2026`; written `MM2_GEO=1 MM2_GEO_CERT=1` in earlier runs), and `mm2`/`stock`
minimap2 (`MM2_GEO` unset or 0). Timings are wall-clock seconds and peak RSS in kB (`/usr/bin/time`); `load=`
is the 1-minute load average when each run started.

## `cleanroom/` — re-run from a fresh clone of this branch

Outputs of `scripts/run_all.sh` run on a fresh clone from GitHub (commit d899b11), with the conda environment
from `environment.yml`, executables built by `build_tools.sh` from the public tags (certified mode as the
default) and inputs downloaded by `fetch_data.sh`; every step ran from an empty environment and all 39 steps
passed (`status.tsv`). Inputs and the work directory were on network storage; chr2–chr10 of script 24 were
re-timed alone (`results/genomewide/summary_retime.txt`), and scripts 10 and 25, whose 32-thread and
disk-output timings are limited by network writes, were re-timed with inputs and output on local disk
(`localdisk/`, the values used in Supplementary Tables S12 and S20).

## `submission/` — first submission

`results_*.txt` written by scripts 01–17 into `$WORKDIR/hifi_bench/` (Table 1, Table 3, Supplementary
Tables S1–S14 and the stress tests); `results_2ndcaller_bcftools.txt`, `results_3rdcaller_dv.txt` (Table S5,
now script 32) and `results_KILLER_f1.txt` (Table S8, now script 33).

## `revision/` — revised manuscript

| folder | content | table | script |
|---|---|---|---|
| `modes_timing/` | `identity.tsv` (certified vs minimap2 hashes), `chr22_timing.tsv` (chr22/chr14, median of 3), `timing.tsv` and `amdahl.tsv` (ONT R9.4.1 genome-wide subset) | Table 2 (end-to-end), S22 | 23 |
| `dp_stage/` | `dpstage.tsv`: wall, total CPU and DP-stage CPU time per data set and mode; `ont_r9_genomewide_amdahl.tsv`: the same for the ONT R9.4.1 genome-wide subset | Table 2 (DP stage), S21 | 26 |
| `gap_rounds/` | band-doubling rounds and joint (d, E) statistics | S15 | 18 |
| `discordant_reads/` | discordant reads by severity and GIAB stratification; `eval_geo.txt`: simulated truth | S16 | 19, 20 |
| `callset_diff/` | `sm_*`/`smhc_*`: small-variant comparison of the heuristic and certified callsets against minimap2 (whole chr22 / GIAB high-confidence); `sv_vcf/`: Sniffles2 SV calls for the nine chr22 BAMs | S17 | 21 |
| `sv_truth_grch37/` | SV calls against GIAB Tier1 v0.6 for minimap2 / heuristic / certified | S17 (bottom) | 22 |
| `genomewide/` | full-coverage HG002 HiFi per chromosome: heuristic series (`wg_perchr.tsv`, `rerun_chr11.tsv`); both modes with the released executables, used for the certified series (`wg_modes_release.tsv`, identical to `cleanroom/results/genomewide/perchr.tsv`); the 10 chr8 reads whose certified output differs (`chr8_discordant_release.paf`: nine equal-score ties, one higher-scoring; `chr8_check.txt`: the same reads in a development build); chr1-only setting (`chr1only_full.tsv`) | S18 | 24 |
| `speed_checks/` | thread counts, output handling, index scope, bioconda and -O3 builds, rebuild from the public tag `submission-2026` (`pubcheck.tsv`, `relcheck.txt`) | S20 | 25 |
| `input_counts/` | reads, bases, mapped and primary rates per data set | S19 | 29 |
| `ont_r9_genomewide/` | ONT R9.4.1 genome-wide subset: read counts and certified identity | S22 | 23 |
| `ont_r10_genomewide/` | ONT R10.4.1 genome-wide subset (10%, seed 42): read counts, certified identity and the 4 differing PAF lines (2 reads, both higher-scoring in the certified mode), timing (3 reps), DP-stage decomposition | S19, S21, S22 | 27 |
| `cert_oracle/` | per-gap oracle of the certified mode (`summary.tsv`) and exhaustive check of the score bound (`bound_exhaustive.tsv`) | S23 | 30, 31 |
| `modes_development/` | timings of the certified mode during development (score certificate alone, with the one-gap shortcut, with the score-only probe) | S22 | 23 |
