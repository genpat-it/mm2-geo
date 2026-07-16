# mm2-geo

**Geometry-guided band-doubling for the minimap2 gap-filling DP** — a small,
runtime-gated patch on top of [minimap2](https://github.com/lh3/minimap2) v2.30
that accelerates the between-anchor dynamic programming at near-identical
variant-calling accuracy.

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
# geo mode:
MM2_GEO=1 ./mm2geo -ax map-hifi ref.fa reads.fq > out.sam
MM2_GEO=1 MM2_GEO_MARGIN=20 ./mm2geo -ax map-hifi ref.fa reads.fq > out.sam
```

- `MM2_GEO=1` — enable geometry-guided band-doubling.
- `MM2_GEO_MARGIN=<int>` — initial band margin (default 20).
- `MM2_GEO_UNCAP=1` — (experimental) let the band exceed minimap2's own bound on
  large-drift gaps.

**Diagnostic flags** (gated; zero measurable cost when unset — a single
never-taken branch — used to reproduce the paper's mechanism/ablation):

- `MM2_GEO_STATS=1` — print logical DP-cell count + per-gap band statistics to stderr.
- `MM2_GEO_FIXED=1` — ablation: fixed band `b=m`, ignoring the geometry.
- `MM2_GEO_NODOUBLE=1` — ablation: single pass, no retry.

## Results (paper)

- **Faster, less memory**: 1.15–1.47× on the alignment stage (1.15–1.35× HiFi/ONT,
  1.38–1.47× CLR) at lower peak memory (up to −43%).
- **Mechanism, measured**: mm2-geo evaluates ~46% of minimap2's gap-filling
  logical DP cells (a 2.2× reduction, near-identical on HiFi and ONT); ~99.9% of
  gaps accept on the first pass at a median band of 20–21.
- **Accuracy preserved**: small-variant F1 within **0.0011** of minimap2 on GIAB
  HG002 (chr22, chr14; HiFi/ONT), three-caller confirmed on HiFi chr22 (clair3,
  bcftools, DeepVariant); within **0.009** on a real *E. coli* ONT spike-in
  benchmark; structural-variant F1 within **0.008** (Sniffles2 + truvari vs GIAB
  Tier1).
- **Complementary to SIMD acceleration**: the same patch applied to
  [mm2-fast](https://github.com/bwa-mem2/mm2-fast) adds a further 1.17–1.23×
  (version-matched; up to ~1.83× vs minimap2 2.30, a figure that spans minimap2
  versions) at near-identical F1.
- Not reducible to lowering minimap2's `-r`: a fixed small band clips indels
  (INDEL F1 drops; on ONT it is even slower), whereas mm2-geo preserves F1.
- Short reads (Illumina, `-ax sr`) use minimap2's ungapped path and are
  unaffected — the method targets long-read gap DP.

## Applying to mm2-fast / mm2-plus

The change touches only the shared gap-filling loop in `align.c`, so it ports
across the minimap2 family:

```bash
cd mm2-fast   # or mm2-plus
patch -p1 < /path/to/geo.patch    # geo.patch is included in this repo
make
MM2_GEO=1 ./minimap2 ...
```

## Honest note

mm2-geo is **not byte-identical** to stock minimap2 (98.3–99.9% CIGAR identity on
the routine benchmarks): minimap2's own gap alignment is approximate
(`KSW_EZ_APPROX_MAX`), and because the accepted band may exclude paths minimap2
would explore, byte-identical output cannot be guaranteed in general. What
mm2-geo provides is **near-identical downstream variant-calling F1**, verified
against GIAB truth across three callers and extended to a real bacterial genome.

## Reproducibility

Scripts, data manifests and raw outputs that regenerate every table and figure of
the paper are in the companion reproducibility repository (see the paper's Code /
Data Availability). A permanent Zenodo snapshot of the code and the
reproducibility material will be archived with a DOI upon acceptance.

## Credit

Built on minimap2 by Heng Li (https://github.com/lh3/minimap2), MIT license.
mm2-geo modifications: Andrea de Ruvo.
