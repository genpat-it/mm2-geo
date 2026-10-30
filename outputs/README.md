# Raw outputs

Summary files as produced by the runs reported in the paper. Large intermediates (BAM, SAM/PAF, read files)
are not included; the scripts in `scripts/` regenerate them. In these files `fast` and `geo` denote the
heuristic mode (`MM2_GEO=1 MM2_GEO_FAST=1`; `MM2_GEO=1` in `submission-2026`), `cert` the certified mode (`MM2_GEO=1`, the default in `revision-2026`; written `MM2_GEO=1 MM2_GEO_CERT=1` in earlier runs), and `mm2`/`stock`
minimap2 (`MM2_GEO` unset or 0). Timings are wall-clock seconds and peak RSS in kB (`/usr/bin/time`); `load=`
is the 1-minute load average when each run started. Memory summaries labelled `GB` by the older scripts (01, 11) are kB/2^20,
i.e. GiB, as printed in the paper; the scripts now print `GiB`.

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
`localdisk/speedchecks_timing.tsv` and `localdisk/25.log` (script 25: timing table and full log). They are kept for
reference: Supplementary Table S12 comes from `release/results/certified/threads.tsv` (script 37) and the
heuristic-mode rows of Table S20 from `revision/speed_checks/`.

The certified-mode identity checks recorded in `cleanroom/` (executables built at commit c042cdc) and in
`revision/ont_r10_genomewide/stats.txt` (an earlier development build, committed on 8 October) predate the release;
the certified-mode identity and every certified-mode timing reported in the paper are those of `release/`, obtained
with the executables of tag `revision-2026`. Heuristic-mode-only timings come from the series stated in each table:
Table S13 from `submission/`, the first rows of Table S14 from `cleanroom/hifi_bench/` (time and memory from `results_WHOLEGENOME.txt`, whose
concordance line counts all reads, the concordance block from the primary-only `results_WG_CONCORD.txt` of script 12), the heuristic-mode rows of Table S20 from
`revision/speed_checks/` (`submission-2026` executables), those of Table S18 from `release/` (`submission-2026`
executable, in the same series as the certified mode). Supplementary Table S18 comes from `release/results/genomewide/perchr.tsv`, one series with
both modes interleaved; the earlier per-mode series in `revision/genomewide/` are superseded.

## `release/` — certified-mode steps with the executables of tag `revision-2026`

The steps that involve the certified mode (00, 13, 14, 18, 20–24, 26–28, 30, 34; `logs/status.tsv`) and the concordance and ablation scripts 05 and 03, and the discordant-read statistics of script 19, were run with the executables built by `build_tools.sh` from tag
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

Scripts 36–38 (tag `revision-2026`; 32 threads unless stated, and one thread in script 38, which is not timed; timed configurations interleaved within each replicate):
`results/table1_modes/` (script 36: `timing.tsv` and `summary.txt`, Table 1 time and memory of minimap2, certified
and heuristic mode; `identity.tsv`, certified identical to minimap2 on all six chr22/chr14 data sets);
`results/certified/` (script 37, one file per section: `chr1only.tsv` Table S18 note, `speed.tsv` Table S20
certified series, `threads.tsv` Table S12, `margin.tsv` Table S9 times and certified identity, `rband.tsv` Table S8
speed-ups, `gwsub.tsv` Table S14 second series, `cells.tsv` Table S2, `ecoli.tsv` certified identity on *E. coli*;
`progress.txt` records that the margin section was re-run so that both modes are timed in the same series, the
first, certified-only run being `margin_certonly_20261030.tsv`). In `margin.tsv` the only `NO` (ONT, m=5) is one
read aligned to a different path of equal score (same AS and NM). `results/examples/` (script 38): PAF of the six
chr22 HiFi reads of the Supplementary "Two examples" in the three configurations and the per-gap log of read
`m64011_190901_095311/1378531/ccs`. `results/ont_genomewide_sam_identity/` (script 39): `summary.tsv`, SAM
(without `@PG`) of minimap2 and of the certified mode on the R9.4.1 and R10.4.1 genome-wide subsets against GRCh38;
identical on both. `results/ecoli_timing/` (script 40): `timing.tsv` and `summary.txt`, *E. coli* ONT time and
memory of minimap2 and both modes (Supplementary text); it replaces `submission/results_ECOLI_TIME_REFMUT.txt`, a run
not scripted in this branch.

## `submission/` — first submission

`results_*.txt` written by scripts 01–17 into `$WORKDIR/hifi_bench/` (Table 1, Table 3, Supplementary
Tables S1–S14 and the stress tests); `results_2ndcaller_bcftools.txt`, `results_3rdcaller_dv.txt` (Table S5,
now script 32) and `results_KILLER_f1.txt` (Table S8, now script 33).

## `revision/` — revised manuscript

| folder | content | table | script |
|---|---|---|---|
| `modes_timing/` | `identity.tsv` (certified vs minimap2 hashes) and `chr22_timing.tsv` (chr22/chr14, median of 3), script 23; `timing.tsv` and `amdahl.tsv` (ONT R9.4.1 genome-wide subset), script 34 | earlier series, superseded by `release/results/table1_modes/` (Table 2 end-to-end, S22) and `release/results/ont_r9_genomewide/` | 23, 34 |
| `dp_stage/` | `dpstage.tsv`: wall, total CPU and DP-stage CPU time per data set and mode (script 26); `ont_r9_genomewide_amdahl.tsv`: the same for the ONT R9.4.1 genome-wide subset (script 34) | earlier series, superseded by `release/results/dpstage/dpstage.tsv` and `release/results/ont_r{9,10}_genomewide/dpstage.tsv` (Table 2 DP stage, S21) | 26, 34 |
| `gap_rounds/` | band-doubling rounds and joint (d, E) statistics of the heuristic mode on chr22 HiFi/ONT/CLR | earlier series; S15, with the genome-wide subset and the certified-mode block, uses `release/results/gap_rounds_dE.txt` | 18 |
| (see `release/`) | discordant reads by severity and GIAB stratification (`release/results/discordant_reads.txt`) and simulated truth (`release/results/simulated_truth.txt`) | S16 (discordant reads); simulated truth in the text of the Supplementary section on the CIGAR-discordant reads | 19, 20 |
| `callset_diff/` | `sm_*`/`smhc_*`: small-variant comparison of the heuristic and certified callsets against minimap2 (whole chr22 / GIAB high-confidence); `sv_vcf/`: Sniffles2 SV calls for the nine chr22 BAMs | S17 | 21 |
| `sv_truth_grch37/` | SV calls against GIAB Tier1 v0.6 for minimap2 / heuristic / certified | S17 (bottom) | 22 |
| `genomewide/` | full-coverage HG002 HiFi per chromosome: heuristic series (`wg_perchr.tsv`, `rerun_chr11.tsv`); both modes with the released executables, an earlier series (`wg_modes_release.tsv`; the clean-room re-run is `cleanroom/results/genomewide/perchr.tsv`); the 9 chr8 reads whose certified output differs, all equal-score ties, are in `release/results/genomewide/discordant_chr8.paf`; chr1-only setting (`chr1only_full.tsv`) | earlier series; S18 uses `release/results/genomewide/perchr.tsv` | 24 |
| `speed_checks/` | thread counts, output handling, index scope, bioconda and -O3 builds, rebuild from the public tag `submission-2026` (`pubcheck.tsv`, `relcheck.txt`) | S20 | 25 |
| `input_counts/` | reads, bases, mapped reads and primary records per data set (`minor3.tsv`; its `certified` rows are certified-mode counts); the heuristic-mode counts cited in the S19 caption are in `cleanroom/results/input_counts.tsv` | S19 | 29 |
| `ont_r9_genomewide/` | ONT R9.4.1 genome-wide subset: input reads and bases, certified identity | S19 (input reads and bases); mapped counts, identity, timing and DP stage from `release/results/ont_r9_genomewide/` (S19, S21, S22) | 34 |
| `ont_r10_genomewide/` | ONT R10.4.1 genome-wide subset (10%, seed 42), earlier development build: read counts, certified identity (2 reads differ, AS +2), timing (3 reps), DP-stage decomposition | S19 (read counts); timing and identity superseded by `release/results/ont_r10_genomewide/` (S21, S22) | 27 |
| `cert_oracle/` | per-gap oracle of an earlier certified build (`summary.tsv`) and exhaustive check of the score bound (`bound_exhaustive.tsv`) | bound check only; the S23 gap counts are in `release/results/cert_oracle/summary.tsv` | 30, 31 |
