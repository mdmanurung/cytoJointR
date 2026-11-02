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

## Model redesign benchmark (2026-11-02)

Source: `cytoJointR-reproducibility/runs/2026-11-02_redesign-summary/` (tables
named below) and `runs/2026-11-02_redesign-decision/README.md`. 545 replicates,
18 A/B/C cells (J = 3 and 4, 24 donors), 30 seeds per cell. C0 is the
independent Dirichlet baseline, C1 is `model = "m1"`, C2 is `model = "m1dm"`
(the default). The pre-registered rule found no calibrated model; the default
is a post hoc choice (`docs/adr/0002-default-model-m1dm.md`).

**Coverage of 90% intervals** (`calibration.tsv`; sample draws, sample target,
global scores). Mean over cells, range in brackets. Generator A has 14 cells,
B and C have 2 each.

| Score | Model | Generator A | Generator B | Generator C |
|---|---|---|---|---|
| RMS | C0 | 0.93 (0.90-0.97) | 0.93 (0.92-0.93) | 0.93 (0.92-0.94) |
| RMS | C1 | 0.91 (0.87-0.99) | 0.91 (0.91-0.91) | 0.87 (0.87-0.88) |
| RMS | C2 | 0.92 (0.88-0.99) | 0.90 (0.90-0.90) | 0.88 (0.88-0.88) |
| CPS | C0 | 0.80 (0.00-1.00) | 0.11 (0.00-0.21) | 0.90 (0.81-0.99) |
| CPS | C1 | 0.78 (0.51-1.00) | 0.30 (0.24-0.35) | 0.69 (0.59-0.79) |
| CPS | C2 | 0.99 (0.98-1.00) | 0.83 (0.81-0.84) | 0.99 (0.99-0.99) |
| CCS | C0 | 0.61 (0.00-0.98) | 0.00 (0.00-0.01) | 0.69 (0.64-0.74) |
| CCS | C1 | 0.99 (0.83-1.00) | 0.33 (0.27-0.39) | 1.00 (1.00-1.00) |
| CCS | C2 | 1.00 (0.98-1.00) | 0.36 (0.34-0.37) | 1.00 (0.99-1.00) |

**Rare-weak detection** (`auc.tsv`; AUC for telling the rare_weak scenario from
null, sample target, J = 3 / J = 4):

| Detector | AUC |
|---|---|
| C1 P(RMS > eps) | 0.989 / 0.987 |
| C2 P(RMS > eps) | 0.983 / 0.980 |
| C0 P(RMS > eps) | 0.805 / 0.799 |
| COMPASS default FS | 0.918 / 0.912 |
| COMPASS nofilter FS | 0.887 / 0.921 |
| polyICSFlow net_sum (unclipped) | 0.810 / 0.802 |

**RMSE gain** (`calibration.tsv`, global RMS, mean over rare_weak and
rare_moderate at J = 3 and 4): C0 6.39e-4, C1 4.42e-4, C2 4.42e-4. Both
hierarchies are 31% below C0.

**False signs** (`false_sign.tsv`; P(effect > eps) > 0.95 while truth <= 0;
generator B has no rows). CCS under C0: up to 0.995 (mean 0.30 in generator A).
CCS under C1: at most 0.002. CCS under C2: 0. RMS false signs are at most 0.058
for every model. CPS false signs are at most 0.006.

**Runtime** (`runtime.tsv`; median seconds per fit, C0 on 1 CPU, C1 and C2 on
4):

| Scenario | C0 | C1 | C2 | COMPASS default | polyICSFlow |
|---|---|---|---|---|---|
| rare_weak_J3, 24 donors | 2.84 | 11.27 | 17.85 | 10.23 | 5.69 |
| runtime_100_J4, 100 donors (5 seeds) | 20.06 | 85.29 | 127.70 | 74.93 | 26.02 |

The slowest C2 fit at 100 donors took 959 s.

## Held-out generator D: PENDING

**Not run yet. Fill this in after D is opened.** Record: seeds and scenarios
(plan: offset 900000, rare_weak and moderate, 30 seeds), RMS and CPS coverage
per model, false signs, and whether `m1dm` still holds. D can only veto: a
default that fails D is reported with the failure, not re-tuned.

## Main limitation: CCS is not calibrated

**The generator-B RMS and CPS collapse is fixed by `m1dm`; CCS is still not
calibrated, and CPS intervals are conservative.** Under generator B the
pilot's M1 gave RMS coverage 0.23 and CPS 0.07. In the redesign benchmark,
`m1dm` gives RMS 0.90 and CPS 0.83 there (`calibration.tsv`, table above).

- **CCS** covers only 0.34 to 0.37 under `m1dm` and 0.27 to 0.39 under `m1`
  in generator B. No model covered it, so CCS is marked experimental: the
  `experimental` column, one message per `cj_scores()` call, and the
  "(experimental)" plot label. The C0 baseline gives CCS false signs up to
  0.995, so dropping the model is no fix.
- **CPS under `m1dm`** covers 0.98 to 1.00 in generators A and C for a nominal
  0.90, so its intervals are wider than needed. In generator B it covers
  0.81 to 0.84. The cause of the over-coverage is unverified; the Dirichlet
  residual may inflate ratio-type scores.
- Generators A, B and C are simulations. Real data may differ from all three.

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
