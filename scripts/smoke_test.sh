#!/bin/bash
# smoke_test.sh — 30-second pre-flight check that YOUR environment can run the pipeline,
# BEFORE you download GB of GIAB data or launch the multi-hour scripts.
#
# It generates a tiny synthetic reference + reads (deterministic, no downloads), then checks:
#   1. the mm2-geo binary ($MM2GEO) runs;
#   2. MM2_GEO=0 (stock) and MM2_GEO=1 (geo) both produce valid, non-empty alignments;
#   3. samtools reads the output;
#   4. the geometric band is actually engaged under MM2_GEO=1 (MM2_GEO_STATS);
#   5. the certified mode (MM2_GEO=1 MM2_GEO_CERT=1) reproduces the stock output exactly.
# It does NOT check accuracy — only that your build + tools + config are wired correctly.
set -uo pipefail
source "$(dirname "$0")/config.sh"
T="$TMPDIR/smoke.$$"; mkdir -p "$T"
pass=1; ok(){ echo "  [OK]  $1"; }; bad(){ echo "  [FAIL] $1"; pass=0; }

echo "== mm2-geo smoke test =="
echo "MM2GEO   = $MM2GEO"
echo "SAMTOOLS = $SAMTOOLS"
[ -x "$MM2GEO" ] || command -v "$MM2GEO" >/dev/null && ok "mm2-geo binary found" || bad "mm2-geo binary not found (set MM2GEO in config.sh)"
command -v "$SAMTOOLS" >/dev/null 2>&1 || [ -x "$SAMTOOLS" ] && ok "samtools found" || bad "samtools not found (activate the conda env / set SAMTOOLS)"

# --- generate a tiny deterministic reference + reads (no external tools/data) ---
python3 - "$T" <<'PY'
import random,sys
T=sys.argv[1]; random.seed(20260716)
ref="".join(random.choice("ACGT") for _ in range(30000))
open(f"{T}/ref.fa","w").write(">smoke_ctg\n"+"\n".join(ref[i:i+70] for i in range(0,len(ref),70))+"\n")
def mut(s,rate=0.01):
    out=[]
    for c in s:
        r=random.random()
        if r<rate*0.6: out.append(random.choice("ACGT"))            # substitution
        elif r<rate*0.8: continue                                   # deletion
        elif r<rate: out.append(c+random.choice("ACGT"))            # insertion
        else: out.append(c)
    return "".join(out)
with open(f"{T}/reads.fq","w") as f:
    for k in range(60):
        L=random.randint(2000,6000); st=random.randint(0,len(ref)-L)
        seq=mut(ref[st:st+L])
        f.write(f"@r{k}\n{seq}\n+\n{'I'*len(seq)}\n")
print("generated 30kb ref + 60 reads")
PY

# --- 1) stock path (MM2_GEO unset) ---
MM2_GEO=0 "$MM2GEO" -ax map-hifi -t4 "$T/ref.fa" "$T/reads.fq" 2>/dev/null > "$T/stock.sam"
n0=$("$SAMTOOLS" view -c -F0x904 "$T/stock.sam" 2>/dev/null || echo 0)
[ "${n0:-0}" -ge 1 ] && ok "MM2_GEO=0 mapped $n0/60 primary reads" || bad "MM2_GEO=0 produced no alignments"

# --- 2) geo path (MM2_GEO=1) ---
MM2_GEO=1 MM2_GEO_MARGIN="${MM2_GEO_MARGIN:-20}" "$MM2GEO" -ax map-hifi -t4 "$T/ref.fa" "$T/reads.fq" 2>/dev/null > "$T/geo.sam"
n1=$("$SAMTOOLS" view -c -F0x904 "$T/geo.sam" 2>/dev/null || echo 0)
[ "${n1:-0}" -ge 1 ] && ok "MM2_GEO=1 mapped $n1/60 primary reads" || bad "MM2_GEO=1 produced no alignments"

# --- 3) geometric band actually engaged ---
stats=$(MM2_GEO=1 MM2_GEO_STATS=1 "$MM2GEO" -ax map-hifi -t1 "$T/ref.fa" "$T/reads.fq" 2>&1 >/dev/null | grep -c "MM2_GEO_STATS\|MM2_GEO_CELLS")
[ "${stats:-0}" -ge 1 ] && ok "geometric band engaged under MM2_GEO=1 (diagnostics printed)" || bad "MM2_GEO_STATS produced no diagnostics"

# --- 4) certified mode == stock (SAM records, @PG excluded) ---
MM2_GEO=1 MM2_GEO_CERT=1 "$MM2GEO" -ax map-hifi -t4 "$T/ref.fa" "$T/reads.fq" 2>/dev/null > "$T/cert.sam"
if cmp -s <(grep -v '^@PG' "$T/stock.sam") <(grep -v '^@PG' "$T/cert.sam"); then ok "certified mode output identical to stock"
else bad "certified mode output differs from stock (is \$MM2GEO built from tag revision-2026?)"; fi

rm -rf "$T"
echo
if [ "$pass" = 1 ]; then echo "SMOKE TEST PASSED — environment is ready; proceed to the numbered scripts."; exit 0
else echo "SMOKE TEST FAILED — fix config.sh / conda env before running the full pipeline."; exit 1; fi
