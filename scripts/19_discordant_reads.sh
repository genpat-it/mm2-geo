#!/bin/bash
# Localisation and severity of the reads whose heuristic-mode alignment differs from minimap2
# (Supplementary Table S16): score loss, truncation, GIAB v3.3 stratifications.
# Reuses the PAF files written by 18_gap_rounds_dE.sh.
set -uo pipefail
source "$(dirname "$0")/config.sh"; cd "$WORKDIR"
H="$(cd "$(dirname "$0")" && pwd)"
[ -s results/gaplog/hifi_geo0.paf ] || { echo "run 18_gap_rounds_dE.sh first"; exit 1; }
python3 "$H/19_discord_stats.py" results/gaplog strat conf_chr22.bed | tee results/discordant_reads.txt
