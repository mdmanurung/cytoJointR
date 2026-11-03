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
3. **Count model in Stan** (only engine): `m1dm` is the default; `m1` is kept
   for ablation. Both are the M1 log-linear pattern model with pairwise
   interactions; `m1dm` adds a Dirichlet-multinomial sample residual (see
   `docs/adr/0002-default-model-m1dm.md`); design = several stimulations per donor sharing a
   control sample, repeated categorical visits, unbalanced donors allowed. One parent population per fit in v0.1.
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
7. **Adoption API** (added 2026-11-03): a user with a cell data frame gates,
   builds a `cj_data` and gets scores in a few calls. `cj_gate()` (pooled
   control gate on a stated transform), `cj_from_cells()` (keeps per-pattern
   intensity sums), `cj_run_compass()`, `cj_run_polyicsflow()`,
   `cj_fluorescence_scores()` and the thin `cj_quickstart()`.
   `cj_fluorescence_scores()` is the descriptive part of the former M2: NICR,
   Poly-NICR, the frequency/intensity split, FIS and experimental
   intensity-weighted FS/PFS, as plug-in estimates from events with a
   cell-level bootstrap. No new Stan model. openCyto, COMPASS and polyICSFlow
   are `Suggests`.

## Deferred

A fluorescence model (M2: truncated-Gaussian conditional intensities with
posterior uncertainty for NICR, Poly-NICR and FIS; the descriptive versions are in
v0.1, item 7); RTMB engine; group terms
and `cj_group_contrast()`; donor-specific induced change across visits;
per-visit interaction change; continuous-time visits; Shiny app; hierarchical
missing-control prediction; CRAN/Bioconductor build layer.

## Pilot benchmark (in the reproducibility folder)

Generator A (exact M1): J=3,4; 2 stimulations; 2 visits; N=10k; donors 24 and
100; truths null, positive margins, negative margins, dependence-only at fixed
margins, visit change. One generator-B (logistic-normal) misspecification.
Smoke 5 seeds, pilot 30 seeds. Report recovery, score bias, 90% interval
coverage, R-hat/divergences; COMPASS FS/PFS on the same data.
