#!/bin/bash
# Build every executable used by the numbered scripts into $BIN (see config.sh), from public sources:
#   mm2geo             tag revision-2026 (heuristic + certified modes), default Makefile (-O2)
#   mm2geo_submission  tag submission-2026 (heuristic mode as first released)
#   mm2geo_instr       revision-2026 + instrumentation.patch (MM2_GEO_GAPLOG, MM2_GEO_TIMEDP; same output)
#   mm2geo_o3          revision-2026 built with bioconda's minimap2 flags (-O3 ...)
#   mm2fast / mm2fastgeo  mm2-fast commit 14fe36c, unpatched / with geo_mm2fast.patch (optional: needs a
#                      x86-64 CPU with AVX2; skipped with a warning if the build fails)
# scripts/smoke_test.sh checks the build; scripts/00_verify_output_identity.sh checks output identity.
set -euo pipefail
source "$(dirname "$0")/config.sh"
HERE="$REPO"
mkdir -p "$BIN/src"; cd "$BIN/src"
log(){ echo "[build $(date +%T)] $*"; }
REPO=https://github.com/genpat-it/mm2-geo
[ -d mm2-geo ] || git clone -q "$REPO" mm2-geo
git -C mm2-geo fetch -q --tags --force origin

build_tag(){ # $1=tag $2=output name $3=extra make args (optional) $4=patch to apply (optional)
  local d="$BIN/src/build_$2"
  [ -x "$BIN/$2" ] && { log "$2 present"; return; }
  rm -rf "$d"; git -C "$BIN/src/mm2-geo" worktree add -f -q "$d" "$1"
  [ -n "${4:-}" ] && git -C "$d" apply "$4"
  log "build $2 ($1${4:+ + $(basename $4)}${3:+, $3})"
  if [ -n "${3:-}" ]; then make -C "$d" -j8 CFLAGS="$3" minimap2 >/dev/null; else make -C "$d" -j8 minimap2 >/dev/null; fi
  cp "$d/minimap2" "$BIN/$2"
}
build_tag revision-2026   mm2geo
build_tag submission-2026 mm2geo_submission
build_tag revision-2026   mm2geo_instr "" "$HERE/instrumentation.patch"
# bioconda minimap2 recipe: CFLAGS="${CFLAGS} -g -Wall -O3 -Wc++-compat" on top of conda-forge's default CFLAGS
build_tag revision-2026   mm2geo_o3 "-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -g -Wall -O3 -Wc++-compat"

# mm2-fast (optional): commit 14fe36c, unpatched and with geo_mm2fast.patch from tag revision-2026.
# safestringlib (bundled by mm2-fast) predates GCC 14, whose defaults turn implicit declarations into errors:
# the missing <stdlib.h> include is added and those diagnostics are relaxed; mm2-fast's own code is untouched.
MM2FAST_COMMIT=14fe36c100f6c2aab224d000f3903ca5909640cd
build_mm2fast(){ # $1=output name $2=patch (optional)
  local d="$BIN/src/build_$1" s arch=avx2
  [ -x "$BIN/$1" ] && { log "$1 present"; return 0; }
  rm -rf "$d"; git clone -q https://github.com/bwa-mem2/mm2-fast "$d" || return 1
  git -C "$d" checkout -q $MM2FAST_COMMIT && git -C "$d" submodule update -q --init --recursive || return 1
  s="$d/ext/TAL/ext/safestringlib"
  sed -i '1i #include <stdlib.h>' "$s/safeclib/safeclib_private.h"
  sed -i '/^CFLAGS=/s/$/ -Wno-error=implicit-function-declaration -Wno-implicit-function-declaration -Wno-implicit-int -fcommon/' "$s/makefile"
  if [ -n "${2:-}" ]; then git -C "$d" apply "$2" || return 1; fi
  grep -qw avx512bw /proc/cpuinfo && arch=avx512
  log "build $1 (mm2-fast ${MM2FAST_COMMIT:0:7}, arch=$arch${2:+ + $(basename $2)})"
  make -C "$d" -j8 arch=$arch PROG=minimap2.mm2fast all > "$d.log" 2>&1 || return 1
  cp "$d/minimap2.mm2fast" "$BIN/$1"
}
git -C "$BIN/src/mm2-geo" show revision-2026:geo_mm2fast.patch > "$BIN/src/geo_mm2fast.patch"
if build_mm2fast mm2fast && build_mm2fast mm2fastgeo "$BIN/src/geo_mm2fast.patch"; then log "mm2-fast built"
else log "WARNING: mm2-fast build failed, see $BIN/src/build_mm2fast*.log (only the mm2-fast rows cannot be reproduced)"; fi

log "versions:"; for b in mm2geo mm2geo_submission mm2geo_instr mm2geo_o3; do printf '  %-18s %s\n' "$b" "$("$BIN/$b" --version)"; done
