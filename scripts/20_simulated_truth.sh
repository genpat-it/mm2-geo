#!/bin/bash
# Discordant reads against a known truth: pbsim3 reads simulated from chr22 (fetch_data.sh generates them,
# seed 2024). Counts how many bases minimap2 and each mm2-geo mode place at the true reference position, for the
# reads whose primary alignment differs (Supplementary Section on discordant reads; Table S16 text).
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"; mkdir -p results/sim
H="$SCRIPTS"
[ -s "$SIM_READS" ] || { echo "simulated reads missing: run fetch_data.sh with pbsim3 on PATH"; exit 1; }
for m in "stock MM2_GEO=0" "heuristic MM2_GEO=1" "certified MM2_GEO=1:MM2_GEO_CERT=1"; do set -- $m
  env ${2//:/ } "$MM2GEO" -cx map-pb -t"$THREADS" chr22_named.fa "$SIM_READS" > results/sim/$1.paf 2>/dev/null; done
{ echo "heuristic vs minimap2:"; python3 "$H/20_truth_eval.py" results/sim/stock.paf results/sim/heuristic.paf "$SIM_MAF"
  echo "certified vs minimap2:"; python3 "$H/20_truth_eval.py" results/sim/stock.paf results/sim/certified.paf "$SIM_MAF"; } | tee results/simulated_truth.txt
