# Implementation status

Against `docs/SCOPE.md` (v0.1). Status: DONE (implemented and tested),
PARTIAL (implemented, with a named gap), NOT_RUN (written, never executed).
Test files are in `tests/testthat/`.

## Package

| Item | Status | Evidence |
|---|---|---|
| Roxygen docs with running examples | DONE | every export has a block; hand-written Rd converted to roxygen |
| testthat suite, slow Stan tests skipped on CRAN | DONE | `skip_if_no_stan()` in `helper-model.R` |
| Two vignettes | DONE | `vignettes/cytoJointR.Rmd`, `vignettes/comparison.Rmd` |
| `R CMD check` 0 errors / 0 warnings | PARTIAL | last run (2026-11-02, host): 0 errors, 1 warning (`qpdf` missing on host), 1 note (dependency check segfaults loading flowWorkspace). One intermittent test error: `test-model-r2-donor.R:35`, missing cmdstanr CSV in a temp dir; the file passes 15/15 on rerun. Cause not established. |
| NEWS.md, pkgdown config | DONE | `NEWS.md`, `_pkgdown.yml` |
| GitHub Actions check | NOT_RUN | `.github/workflows/R-CMD-check.yaml` exists; the repository has no remote, so CI has never run |
| Stan models in `inst/stan/`, compiled by cmdstanr | DONE | `inst/stan/m1.stan`, `test-model-*.R` |

## In v0.1

| # | Item | Status | Evidence |
|---|---|---|---|
| 1 | `cj_from_counts()` | DONE | `test-input.R`, `test-classes.R` |
| 1 | Events to counts (`cj_from_events()`) | DONE | `test-input.R` ("events, polyICSFlow and the pattern tabulator agree with counts"), `test-sufficient-statistics.R` |
| 1 | `cj_from_gatingset()` | NOT_RUN | test in `test-input.R` skips; flowWorkspace segfaults on the development host even after rebuilding its dependency chain from source |
| 1 | polyICSFlow `extractCellCounts()` input | PARTIAL | `cj_from_polyicsflow()` tested only on a hand-built table (`test-input.R`), not on real polyICSFlow output |
| 2 | Descriptive layer: counts, signed net, display-only clip | DONE | `test-describe.R` |
| 2 | Pattern heatmap, exactly-n and at-least-n bars, 2D scatters | DONE | `test-describe.R` (plots build as ggplot objects); scatter has an example only |
| 3 | M1 count model, fixed effects | DONE | `test-model-r1-fixed.R`, `test-model-design.R`, `test-model-prior.R` |
| 3 | Several stimulations sharing a control; unbalanced donors; missing controls | DONE | `test-model-r3-stims.R` |
| 3 | Donor effects `b_i`, `r_im`, Sigma | DONE | `test-model-r2-donor.R` |
| 3 | Visits, visit x stimulation, donor x visit deviation and its ablation | DONE | `test-model-r4-visits.R` |
| 4 | Scores from aligned draws (RMS, poly, CPS, CCS, marker magnitude, patterns) | DONE | `test-scores.R`, `test-score-engine.R`, `test-signed-attribution.R`, `test-counterfactual.R` |
| 4 | Three-way epsilon probabilities, unavailable contrasts | DONE | `test-scores.R` |
| 4 | Visit change per donor and at equal donor weights | DONE | `test-visit-change.R` |
| 4 | Per-donor score draws exported (`cj_score_draws()`) | DONE | `test-scores.R`, `helper-scores-draws.R` |
| 5 | COMPASS adapter: FS, PFS, `Response()`, filters, heatmap | DONE | `test-compass.R` |
| 5 | COMPASS posterior differences | PARTIAL | not wrapped; `COMPASS::PosteriorDiff()` runs on `cj_compass()$fits` and is not tested here |
| 6 | Static ggplot2 plots | DONE | forest, response heatmap, contribution bars in `test-scores.R` |
| 6 | ComplexHeatmap output | PARTIAL | plots are ggplot2; ComplexHeatmap is only a Suggests entry and is not used |
| - | Pilot benchmark | PARTIAL | 13 of 21 scenarios x 30 seeds run; see "Pilot results" below |

## Known blockers

* `cj_from_gatingset()` has never run. flowWorkspace segfaults on the
  development host, even after its dependency chain was rebuilt from source.
* The polyICSFlow adapter is tested only on a hand-built table.
* CI has never run because the repository has no remote.
* The default priors in `cj_default_priors()` are illustrative smoke-test
  values, not validated.

## Pilot results (2026-11-02)

Source: `cytoJointR-reproducibility/runs/2026-11-02_pilot/` (its README pins
the commits and deviations). 30 seeds per scenario; coverage of 90% intervals,
binomial Monte Carlo SE about 0.05 at 30 seeds. One pilot, not confirmation.

* **Generator A (data from M1 itself), RMS coverage 0.90-0.99** across null,
  positive margins, negative margins, dependence-only and visit change, at 24
  and 100 donors, J = 3. Bias within about 3e-3.
* **Under-coverage to watch:** CPS 0.82 and CCS 0.80 (SE 0.07) for
  dependence-only at 24 donors; CCS 0.84 for null. Observation, not conclusion.
* **False sign rate** (P(delta > epsilon) > 0.95 while truth <= 0): 0 to 0.007
  for generator A.
* **Donor-by-visit ablation:** dropping `u_iv` collapses visit-change coverage
  (donor-level RMS 0.03 vs 0.89-0.90 for the full model).
* **COMPASS on the same data (seeds 1-5):** FS/PFS saturate near 1 for any
  positive effect; rank correlation with true RMS 0.56-0.79 (-0.14 at J = 4).
  cytoJointR RMS correlates 0.94-1.00 with truth.
* **Mixing:** 0-5 of 30 replicates per scenario flagged (R-hat > 1.01, ESS < 400
  or a divergence), all kept in the metrics. Dependence-only J = 4 hit max
  treedepth in 29 of 30 fits; treedepth was not a flag criterion.
* **Fit time:** median 30-150 s at 24 donors, 500-1300 s at 100 donors (N = 10 000).
* **Not run:** J = 4 at 100 donors; J = 4 at 24 donors for null, negative
  margins and visit change.

## Main limitation: donor heterogeneity beyond the model

**Under generator B (logistic-normal pattern probabilities, heterogeneous
donors), coverage collapses: RMS 0.23 (SE 0.08), CPS 0.07, CCS 0.14, and the
CCS false sign rate is 0.19.** M1 gives donors one shared basal and one induced
effect with unit loadings across markers, so donor variation that is
marker-specific or pattern-specific is pooled away and intervals are far too
narrow. This is the over-pooling risk the design bundle named (doc 14). Real
data will usually look more like generator B than generator A. Candidate
fixes, not yet built: marker-specific donor effects, or sample-level
overdispersion (for example a per-sample pattern random effect).

## Parameterisation

Default is noncentered. `centered = TRUE` is a speed option: it funnels when a
donor variance is near zero (pilot: R-hat up to 1.4 for null truths; 1.23 when
`donor_visit = TRUE` meets data without donor-by-visit variation). In the pilot
noncentered fits at 100 donors (median 503-643 s) were not slower than centered
ones (1020-1322 s, different scenarios).

## Deferred (from SCOPE)

M2 conditional fluorescence (NICR, Poly-NICR, FIS); RTMB engine; group terms
and `cj_group_contrast()`; donor-specific induced change across visits;
per-visit interaction change; continuous-time visits; Shiny app; hierarchical
missing-control prediction; CRAN/Bioconductor build layer.
