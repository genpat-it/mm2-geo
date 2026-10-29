#!/bin/bash
# Reproduce every table of the paper in one go: build the executables, download the data, run the smoke test
# and then every numbered script in dependency order, one log per step in $WORKDIR/logs/ and a summary
# (exit code and run time of each step) in $WORKDIR/logs/status.tsv. Run it from the activated conda
# environment (`conda activate mm2-geo-repro`) with WORKDIR set.
#
#   bash scripts/run_all.sh                 # everything, sequentially (~2 days; ~400 GB)
#   bash scripts/run_all.sh --core          # without the genome-wide scripts 24, 27 and 34 (~1 day; ~150 GB)
#   bash scripts/run_all.sh --parallel      # genome-wide branch in parallel with the core scripts (~1.5 days)
#
# Already completed steps are skipped on a re-run (delete their line from status.tsv to repeat one).
# Timings assume an otherwise idle machine; --parallel shares the CPUs between two 32-thread runs, which
# leaves the speed-up ratios essentially unchanged but makes absolute times longer.
set -uo pipefail
SCRIPTS="$(cd "$(dirname "$0")" && pwd)"; source "$SCRIPTS/config.sh"
CORE_ONLY=0; PARALLEL=0
for a in "$@"; do case $a in --core) CORE_ONLY=1;; --parallel) PARALLEL=1;;
  *) echo "usage: bash scripts/run_all.sh [--core] [--parallel]"; exit 2;; esac; done
[ -d "$WORKDIR" ] || mkdir -p "$WORKDIR"
L="$WORKDIR/logs"; mkdir -p "$L"; S="$L/status.tsv"; touch "$S"

done_ok(){ grep -q "^$1	exit=0" "$S"; }
step(){ # name command...: run unless already completed, log, record exit code and minutes
  local n=$1; shift; done_ok "$n" && { echo "[skip] $n (already done)"; return 0; }
  echo "[run ] $n  (log: $L/$n.log)"; local t0=$(date +%s)
  ( "$@" ) > "$L/$n.log" 2>&1; local rc=$?
  echo -e "$n\texit=$rc\t$(( ($(date +%s)-t0)/60 )) min\t$(date '+%F %T')" >> "$S"
  [ $rc = 0 ] || echo "[FAIL] $n exit=$rc (see $L/$n.log)"; return $rc; }
num(){ local f; f=$(ls "$SCRIPTS"/${1}_* | grep -v '_stratify\|_gap_stats\|_discord_stats\|_truth_eval' | head -1)
  case $f in */31_*.py) step "$1" python3 "$f" "$WORKDIR/results/cert_oracle/bound_exhaustive.tsv";;  # table also to a file
    *.py) step "$1" python3 "$f";; *) step "$1" bash "$f";; esac; }
check15(){ # regenerate the E. coli spike-in (script 15, fixed seed) and compare it with the shipped artifacts
  local d="$WORKDIR/check15"; rm -rf "$d"; mkdir -p "$d/ecoli_real2/flye" "$d/ecoli_vc"
  cp "$WORKDIR/ecoli_real2/flye/assembly.fasta" "$d/ecoli_real2/flye/" && cd "$d" &&
  python3 "$SCRIPTS/15_ecoli_plant_variants.py" &&
  cmp ecoli_vc/ref_mut.fa "$REPO/data/ecoli_ref_mut.fa" &&
  diff <(grep -v '^#' ecoli_vc/truth.vcf) <(zcat "$REPO/data/ecoli_spikein_truth.vcf.gz" | grep -v '^#') &&
  echo "script 15 output IDENTICAL to data/"; }

step build_tools bash "$SCRIPTS/build_tools.sh" || exit 1
step fetch_core  bash "$SCRIPTS/fetch_data.sh" core || exit 1
step smoke_test  bash "$SCRIPTS/smoke_test.sh" || exit 1

genomewide(){ step fetch_genomewide bash "$SCRIPTS/fetch_data.sh" genomewide && { num 27; num 34; num 24; }; }
core(){
  for n in 00 29 18 19 20 23 26 02 03 05 06; do num $n; done
  step 15 check15
  for n in 16 21 22 07 01 04 09 17 08 10 11 12 13 14 25 28 30 31 32 33; do num $n; done; }

if [ $CORE_ONLY = 1 ]; then core
elif [ $PARALLEL = 1 ]; then genomewide & core; wait
else core; genomewide; fi

echo; echo "== summary ($S) =="; cat "$S"
fails=$(awk -F'\t' '$2!="exit=0"' "$S" | wc -l)
[ "$fails" = 0 ] && echo "ALL STEPS PASSED — results in $WORKDIR/hifi_bench and $WORKDIR/results" || echo "$fails step(s) failed"
