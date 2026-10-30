# mm2-geo

**Geometry-guided band sizing for the minimap2 gap-filling DP** — a small,
runtime-gated patch on top of [minimap2](https://github.com/lh3/minimap2) v2.30
that reduces the between-anchor dynamic programming. Its default **certified mode** preserves
minimap2's gap scores (output identical to minimap2's except for rare equal-score ties); an opt-in **heuristic mode** trades exact equivalence for speed on some data.

**How the algorithm works (illustrated, with an interactive bound calculator):
<https://genpat-it.github.io/mm2-geo/>**

> mm2-geo is minimap2 v2.30 (Li, 2018; MIT license, preserved in `LICENSE.txt`)
> plus a compact, self-contained change in `align.c`. With `MM2_GEO` unset the executable produces
> **alignment records byte-identical to stock minimap2** (the binary carries the
> gated code, so it is *output*-identical, not a byte-identical executable);
> verified by a record-level `diff` after normalising SAM headers. Original
> minimap2 README: see the upstream repo.

## Idea

Between two chained anchors, minimap2 fills the gap with a banded DP whose width
is a single global heuristic. mm2-geo instead sizes the band **per gap** from the
local geometry `band = |Δquery − Δref| + margin` and **doubles it adaptively**
under a *boundary-contact test* (`max|ΣI − ΣD| < band`) — a heuristic
band-sufficiency signal, **not** a proof of global optimality — never exceeding
minimap2's own band. Clean gaps get a tiny band (fast); the rare hard gap widens
up to minimap2's cap. By Lemma 1 (paper), the geometric floor `|Δq − Δr|` is what
makes the acceptance test operate in its valid regime.

## Build & use

```bash
make
cp minimap2 mm2geo   # the Makefile emits ./minimap2; run it as ./mm2geo below
# stock minimap2 behaviour (output-identical to minimap2 2.30):
./mm2geo -ax map-hifi ref.fa reads.fq > out.sam
# certified mode (recommended; preserves minimap2's gap scores, output identical except rare equal-score ties):
MM2_GEO=1 ./mm2geo -ax map-hifi ref.fa reads.fq > out.sam
# heuristic mode (opt-in; faster on some data, not output-identical):
MM2_GEO=1 MM2_GEO_FAST=1 ./mm2geo -ax map-hifi ref.fa reads.fq > out.sam
```

- `MM2_GEO=1` — certified mode (the default): the band is sized from the anchor
  geometry and accepted only when the in-band score is provably at least that of any
  path leaving it; otherwise the gap is realigned with the smallest band the
  bound certifies and tested again. Includes an exact one-gap shortcut and a score-only first pass on
  gaps of at least 256 bp (`MM2_GEO_CERT_SCOREPROBE=<int>`, 0 disables). Output is
  identical to minimap2's, except for rare equal-score ties (see below).
- `MM2_GEO_FAST=1` — heuristic mode: geometry-sized band, doubled when the returned
  path touches the band boundary (a heuristic acceptance test). This is the mode of
  the first release (tag `submission-2026`, where it was selected by `MM2_GEO=1` alone).
- `MM2_GEO_CERT=1` — accepted for compatibility with earlier commands; no effect.
- `MM2_GEO_ONEGAP=1` — with `MM2_GEO_FAST=1`, adds the exact one-gap shortcut alone.
- `MM2_GEO_MARGIN=<int>` — initial band margin (default 20).

**Diagnostic flags** (gated; zero measurable cost when unset — a single
never-taken branch — used to reproduce the paper's mechanism/ablation):

- `MM2_GEO_STATS=1` — print logical DP-cell count + per-gap band statistics to stderr.
- `MM2_GEO_FIXED=1` — ablation: fixed band `b=m`, ignoring the geometry.
- `MM2_GEO_NODOUBLE=1` — ablation: single pass, no retry.

## Results (paper)

Same executable with only environment variables changing, 32 threads, chr22 GIAB HG002 unless stated (median of three runs; genome-wide: all 11.7 million HG002 HiFi reads, both modes in one series):

| data | heuristic mode | certified mode (preserves gap scores) |
|---|---|---|
| HiFi | 1.39× (1.32× genome-wide, 53×) | 1.25× (1.33× genome-wide, 53×) |
| ONT (chr22, chr14) | 1.16–1.19× | 1.07–1.19× |
| CLR | 1.41× | 1.11× |

- The DP stage itself becomes 1.6–2.1× faster in the heuristic mode and 1.16–2.6×
  faster in the certified mode; the rest of the run is unchanged, so the
  end-to-end gain follows the share of run time spent in DP (Amdahl's law).
- The heuristic mode keeps small-variant F1 unchanged in GIAB high-confidence regions
  but can accept lower-scoring alignments in repetitive sequence; use the
  certified mode for structural-variant or repeat analyses and validated pipelines.
- Complementary to SIMD acceleration: a port to
  [mm2-fast](https://github.com/bwa-mem2/mm2-fast) is provided (`geo_mm2fast.patch`, for the
  genomic presets; spliced alignment is not covered):
  on chr22 its certified mode is byte-identical to mm2-fast and 1.17× (HiFi), 1.08× (ONT)
  and 1.08× (CLR) faster than it; the heuristic mode 1.23×, 1.18× and 1.33×.
- Short reads (Illumina, `-ax sr`) keep minimap2's own gap filling and are
  unaffected.

## Applying to minimap2 and mm2-fast

`geo.patch` applies to stock minimap2 v2.30 (`align.c` only). It is the minimal production patch:
heuristic and certified modes without the diagnostic flags (`MM2_GEO_STATS`, `MM2_GEO_FIXED`,
`MM2_GEO_NODOUBLE`), which are in this tree (together with a DP-cell counter in `ksw2_extd2_sse.c`).
mm2-fast restructures the
gap-filling loop, so it needs its own port, `geo_mm2fast.patch` (heuristic and certified modes),
which applies to mm2-fast commit `14fe36c` and covers the genomic presets (map-hifi, map-ont,
map-pb, ...; not spliced alignment):

```bash
git clone https://github.com/bwa-mem2/mm2-fast && cd mm2-fast
git checkout 14fe36c100f6c2aab224d000f3903ca5909640cd
git apply /path/to/geo_mm2fast.patch
make            # see the mm2-fast README for its build requirements
MM2_GEO=1 ./minimap2 ...                  # certified mode (output identical to mm2-fast)
MM2_GEO=1 MM2_GEO_FAST=1 ./minimap2 ...   # heuristic mode
```

## Honest note

The **heuristic mode** is not byte-identical to stock minimap2 (98.3–99.9% CIGAR
identity): its boundary-contact test can accept a lower-scoring path when a
band-leaving insertion–deletion detour is replaced by mismatches, which happens
mostly in repetitive sequence. The **certified mode** replaces that test with a
score bound and preserves minimap2's gap scores; its guarantee is relative to
minimap2, whose own band, z-drop and chaining remain heuristics. Each accepted band carries
a proof that no wider band (up to minimap2's) scores higher, and the returned traceback is
checked to attain that score (otherwise the gap is realigned with minimap2's band), so
the certified gap score equals minimap2's. The only possible read-level differences are
equal-score alternative tracebacks, which can also change minimap2's CIGAR-based z-drop
decision downstream.

Versions: tag `submission-2026` = heuristic mode only (selected by `MM2_GEO=1`); tag `revision-2026` =
certified mode as the default (`MM2_GEO=1`) and the unchanged heuristic mode (`MM2_GEO=1 MM2_GEO_FAST=1`).

## Reproducibility

Scripts, data manifests and raw outputs that regenerate every table and figure of
the paper are in the companion reproducibility repository (see the paper's Code /
Data Availability). A permanent Zenodo snapshot of the code and the
reproducibility material will be archived with a DOI upon acceptance.

## Credit

Built on minimap2 by Heng Li (https://github.com/lh3/minimap2), MIT license.
mm2-geo modifications: Andrea de Ruvo.
