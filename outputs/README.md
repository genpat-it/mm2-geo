# Raw outputs

Summary files as produced by the runs reported in the paper. Large intermediates (BAM, SAM/PAF, read files)
are not included; the scripts in `scripts/` regenerate them. In these files `fast` and `geo` denote the
heuristic mode (`MM2_GEO=1 MM2_GEO_FAST=1`; `MM2_GEO=1` in `submission-2026`), `cert` the certified mode (`MM2_GEO=1`, the default in `revision-2026`; written `MM2_GEO=1 MM2_GEO_CERT=1` in earlier runs), and `mm2`/`stock`
minimap2 (`MM2_GEO` unset or 0). Timings are wall-clock seconds and peak RSS in kB (`/usr/bin/time`); `load=`
is the 1-minute load average when each run started.

## `cleanroom/` — re-run from a fresh clone of this branch

Outputs of `scripts/run_all.sh` run on a fresh clone from GitHub, with the conda environment from
`environment.yml`, executables built by `build_tools.sh` from the public tags (certified mode as the default)
and inputs downloaded by `fetch_data.sh`. The clone and the build were made at commit c042cdc on 2026-10-23;
the script fixes of commits e9aa046 (heuristic mode via `MM2_GEO_FAST=1` also where `MM2_GEO` is set from a
loop variable) and d899b11 (`fetch_data.sh`: the R9.4.1 stream ends on SIGPIPE by design) were pulled into the
clone and picked up when the run was restarted, so the outputs correspond to the scripts at d899b11. Every
step ran from an empty environment and all 39 steps passed (`status.tsv`). Inputs and the work directory were
on network storage. `results/genomewide/summary_retime.txt` is the manual re-timing of partitions chr2–chr10
of script 24, run alone after their first timing had overlapped with other jobs (`perchr.tsv` holds the
re-timed rows). Scripts 10 and 25, whose 32-thread and disk-output timings are limited by network writes,
were re-timed with inputs and output on local disk: `localdisk/results_THREADSCALING.txt` (script 10),
`localdisk/speedchecks_timing.tsv` and `localdisk/25.log` (script 25: timing table and full log); these are
the values used in Supplementary Tables S12 and S20.

The certified-mode identity checks recorded in `cleanroom/` (the chr8 partition of script 24 and script 27) and in
`revision/ont_r10_genomewide/stats.txt` were made with executables built at commit c042cdc; the certified-mode
identity reported in the paper is that of `release/`, obtained with the executables of tag `revision-2026`.
The timings of Supplementary Table S18 come from `revision/genomewide/wg_perchr.tsv` with `rerun_chr11.tsv`
(heuristic-mode series) and `revision/genomewide/wg_modes_release.tsv` (certified-mode series);
`release/results/genomewide/perchr.tsv` is a later run, used for the per-read identity only.

## `release/` — certified-mode steps with the executables of tag `revision-2026`

The steps that involve the certified mode (00, 13, 14, 18, 20–24, 26–28, 30, 34; `logs/status.tsv`) and the
concordance and ablation scripts 05 and 03 were run with the executables built by `build_tools.sh` from tag
`revision-2026`; these are the outputs behind every certified-mode number in the revised manuscript, and the
concordance and CIGAR-ablation values of Table 3 and Supplementary Table S6 (`hifi_bench/`). Highlights:
`results/identity/summary.txt` and `results/modes/identity.tsv` (certified identical to minimap2 on chr22
HiFi/ONT/CLR, chr14 ONT and the genome-wide HiFi subset, `-cx` and `-ax`; for the subset the run was a byte
comparison, `results/gwsubset_identity/summary.txt`, since added to script 23), `results/genomewide/perchr.tsv` (full-coverage HiFi: identical on every partition except chr8, whose
9 differing reads, all equal-score ties, are in `discordant_chr8.paf`), `results/ont_r9_genomewide/` and
`results/ont_r10_genomewide/` (identical). `results/hifi53x_alignment_rate/summary.txt`: alignment rate of stock minimap2 on all
11,678,114 reads of the full-coverage HiFi data (Supplementary Table S19; script 35, run with the stock minimap2 v2.30
executable, whose output equals that of mm2-geo with `MM2_GEO` unset). Large intermediates (BAM, PAF, Clair3 working directories, truvari
VCFs) are not included.

## `submission/` — first submission

`results_*.txt` written by scripts 01–17 into `$WORKDIR/hifi_bench/` (Table 1, Table 3, Supplementary
Tables S1–S14 and the stress tests); `results_2ndcaller_bcftools.txt`, `results_3rdcaller_dv.txt` (Table S5,
now script 32) and `results_KILLER_f1.txt` (Table S8, now script 33).

## `revision/` — revised manuscript

| folder | content | table | script |
|---|---|---|---|
| `modes_timing/` | `identity.tsv` (certified vs minimap2 hashes) and `chr22_timing.tsv` (chr22/chr14, median of 3), script 23; `timing.tsv` and `amdahl.tsv` (ONT R9.4.1 genome-wide subset), script 34 | Table 2 (end-to-end), S22 | 23, 34 |
| `dp_stage/` | `dpstage.tsv`: wall, total CPU and DP-stage CPU time per data set and mode (script 26); `ont_r9_genomewide_amdahl.tsv`: the same for the ONT R9.4.1 genome-wide subset (script 34) | Table 2 (DP stage), S21 | 26, 34 |
| `gap_rounds/` | band-doubling rounds and joint (d, E) statistics | S15 | 18 |
| (see `release/`) | discordant reads by severity and GIAB stratification (`release/results/discordant_reads.txt`) and simulated truth (`release/results/simulated_truth.txt`) | S16 | 19, 20 |
| `callset_diff/` | `sm_*`/`smhc_*`: small-variant comparison of the heuristic and certified callsets against minimap2 (whole chr22 / GIAB high-confidence); `sv_vcf/`: Sniffles2 SV calls for the nine chr22 BAMs | S17 | 21 |
| `sv_truth_grch37/` | SV calls against GIAB Tier1 v0.6 for minimap2 / heuristic / certified | S17 (bottom) | 22 |
| `genomewide/` | full-coverage HG002 HiFi per chromosome: heuristic series (`wg_perchr.tsv`, `rerun_chr11.tsv`); both modes with the released executables, used for the certified series (`wg_modes_release.tsv`; the clean-room re-run is `cleanroom/results/genomewide/perchr.tsv`); the 9 chr8 reads whose certified output differs, all equal-score ties, are in `release/results/genomewide/discordant_chr8.paf`; chr1-only setting (`chr1only_full.tsv`) | S18 | 24 |
| `speed_checks/` | thread counts, output handling, index scope, bioconda and -O3 builds, rebuild from the public tag `submission-2026` (`pubcheck.tsv`, `relcheck.txt`) | S20 | 25 |
| `input_counts/` | reads, bases, mapped reads and primary records per data set (`minor3.tsv`; its `certified` rows are certified-mode counts) | S19 | 29 |
| `ont_r9_genomewide/` | ONT R9.4.1 genome-wide subset: read counts and certified identity | S19, S21, S22 | 34 |
| `ont_r10_genomewide/` | ONT R10.4.1 genome-wide subset (10%, seed 42): read counts, certified identity (identical to minimap2), timing (3 reps), DP-stage decomposition | S19, S21, S22 | 27 |
| `cert_oracle/` | per-gap oracle of the certified mode (`summary.tsv`) and exhaustive check of the score bound (`bound_exhaustive.tsv`) | S23 | 30, 31 |
