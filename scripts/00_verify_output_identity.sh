#!/bin/bash
set -e
source "$(dirname "$0")/config.sh"
cd "$MM2GEO_SRC"
Q=/tmp/hifi_sub.fq; IDX=/tmp/c22.hifi.mmi
NEW=./mm2geo                 # current (instrumented) binary, already built
BASE=../minimap2_src/minimap2
# build OLD geo (pre-instrumentation) from committed git version -> "$MM2GEO"
git stash push -m instr align.c ksw2_extd2_sse.c >/tmp/stash.log 2>&1
make >/tmp/build_old.log 2>&1
cp -f minimap2 "$MM2GEO"
git stash pop >>/tmp/stash.log 2>&1
make >/tmp/build_new.log 2>&1     # restore instrumented objects
cp -f minimap2 mm2geo
OLD="$MM2GEO"
sam(){ grep -v '^@' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' | sort; }
{
echo "== A) new MM2_GEO=0  vs  baseline minimap2 (OFF gating) =="
diff <(MM2_GEO=0 $NEW -ax map-hifi -t8 $IDX $Q 2>/dev/null | sam) \
     <($BASE       -ax map-hifi -t8 $IDX $Q 2>/dev/null | sam) >/dev/null && echo "IDENTICAL ✓" || echo "DIFFER ✗"
echo "== B) new MM2_GEO=1  vs  OLD geo MM2_GEO=1 (geo path unchanged) =="
diff <(MM2_GEO=1 $NEW -ax map-hifi -t8 $IDX $Q 2>/dev/null | sam) \
     <(MM2_GEO=1 $OLD -ax map-hifi -t8 $IDX $Q 2>/dev/null | sam) >/dev/null && echo "IDENTICAL ✓" || echo "DIFFER ✗"
echo "== C) new MM2_GEO=1  vs  new MM2_GEO=0  (geo actually changes something) =="
diff <(MM2_GEO=1 $NEW -ax map-hifi -t8 $IDX $Q 2>/dev/null | sam) \
     <(MM2_GEO=0 $NEW -ax map-hifi -t8 $IDX $Q 2>/dev/null | sam) >/dev/null && echo "identical (unexpected)" || echo "differ (expected — geo IS active) ✓"
echo "### VERIFY-DONE"
} > /tmp/verify.out 2>&1
