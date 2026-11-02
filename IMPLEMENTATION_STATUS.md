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
| `R CMD check` 0 errors / 0 warnings | DONE | see the check result in the commit log and report |
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
| - | Pilot benchmark | NOT_RUN | running separately; no result is claimed here |

## Known blockers

* `cj_from_gatingset()` has never run. flowWorkspace segfaults on the
  development host, even after its dependency chain was rebuilt from source.
* The polyICSFlow adapter is tested only on a hand-built table.
* The pilot benchmark runs in the reproducibility folder. Parameter recovery,
  score bias and interval coverage at scale are not claimed.
* CI has never run because the repository has no remote.
* The default priors in `cj_default_priors()` are illustrative smoke-test
  values, not validated.

## Deferred (from SCOPE)

M2 conditional fluorescence (NICR, Poly-NICR, FIS); RTMB engine; group terms
and `cj_group_contrast()`; donor-specific induced change across visits;
per-visit interaction change; continuous-time visits; Shiny app; hierarchical
missing-control prediction; CRAN/Bioconductor build layer.
