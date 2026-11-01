# cytoJointR v0.1 scope

Decided 2026-11-01 in a grilling session over the unified design bundle
(`cytojointR_unified_implementation_bundle_2026-11-01.zip`, which supersedes
the earlier bundle). Vocabulary: `GLOSSARY.md`. COMPASS decision: `docs/adr/0001`.
Where this file and the bundle disagree, this file wins for v0.1.

## Package

- Name `cytoJointR`; maintainer Mikhael Manurung; MIT licence.
- Full R package: roxygen2 docs with running examples, testthat suite (slow
  Stan tests skipped on CRAN), two vignettes (count workflow; comparison with
  COMPASS and polyICSFlow), `R CMD check` 0 errors / 0 warnings, NEWS.md,
  pkgdown config, GitHub Actions check.
- Distribution: GitHub / r-universe. Stan models in `inst/stan/`, compiled by
  cmdstanr on first use; cmdstanr via `Additional_repositories`.
- Runtime: conda env `R4_51` (R 4.5.1, CmdStan 2.39) plus a project library in
  the reproducibility folder for COMPASS and a working flowWorkspace.

## In v0.1

1. **Input layer**: `cj_from_counts()` (exact integer pattern counts, parent N,
   pairing table), events-to-counts via the existing compressor,
   `cj_from_gatingset()`, and acceptance of polyICSFlow `extractCellCounts()`
   output. Full validation per bundle doc 01: control samples ingested once,
   all 2^J patterns, no silent drops.
2. **Descriptive layer** (polyICSFlow parity, own implementation): pattern
   counts and percentages, signed background subtraction, clipped values only
   as a display-only column, pattern heatmap, exactly-n / at-least-n bars,
   per-pattern 2D scatters.
3. **M1 count model in Stan** (only engine): log-linear pattern model with
   pairwise interactions; design = several stimulations per donor sharing a
   control sample, repeated categorical visits, unbalanced donors allowed.
   - Fixed: marker main effects; per-stimulation induced effects; visit and
     visit x stimulation effects; pairwise interactions theta0 and per-
     stimulation delta (shared across visits).
   - Donor: basal b_i and per-stimulation induced r_im shared across visits,
     Sigma (1+M)x(1+M); plus donor x visit basal deviation. Ablation: without it.
   - Missing control sample at a visit: that visit's net scores `unavailable`;
     stimulated data still enters the fit.
4. **Scores from aligned draws**: RMS, absolute poly, CPS, poly_marginal, CCS,
   marker magnitude, exact-pattern change; three-way epsilon probabilities;
   visit change per donor and at equal donor weights; per-donor score draws
   exported for downstream outcome analysis.
5. **COMPASS adapter** `cj_compass()`: native FS, PFS, `Response()`,
   posterior differences, mean-gamma heatmap, with COMPASS's own default
   filters (reported). COMPASS in `Suggests`.
6. Heatmaps and plots as static ggplot2 / ComplexHeatmap.

## Deferred

M2 conditional fluorescence (NICR, Poly-NICR, FIS); RTMB engine; group terms
and `cj_group_contrast()`; donor-specific induced change across visits;
per-visit interaction change; continuous-time visits; Shiny app; hierarchical
missing-control prediction; CRAN/Bioconductor build layer.

## Pilot benchmark (in the reproducibility folder)

Generator A (exact M1): J=3,4; 2 stimulations; 2 visits; N=10k; donors 24 and
100; truths null, positive margins, negative margins, dependence-only at fixed
margins, visit change. One generator-B (logistic-normal) misspecification.
Smoke 5 seeds, pilot 30 seeds. Report recovery, score bias, 90% interval
coverage, R-hat/divergences; COMPASS FS/PFS on the same data.
