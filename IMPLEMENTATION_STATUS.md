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
| `R CMD check` 0 errors / 0 warnings | PARTIAL | last run (2026-11-03, host, branch `feat/adoption-api`): 0 errors, 1 warning (`qpdf` missing on host), 2 notes (cannot verify current time; dependency check segfaults loading flowWorkspace). Tests 1583 pass, 3 skip. One intermittent test error: `test-model-r2-donor.R:35`, missing cmdstanr CSV in a temp dir; the file passes 15/15 on rerun. Cause not established. |
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
| 4 | CCS support gate (`ccs_gate`): `not_estimable` with a reason below 20 polyfunctional cells or posterior CV 0.5 of the control's polyfunctional mass | DONE | `test-ccs-gate.R`; thresholds from `scripts/55_calibrate-ccs-gate.R` (`runs/2026-11-03_ccs-gate/`) |
| 4 | Three-way epsilon probabilities, unavailable contrasts | DONE | `test-scores.R` |
| 4 | Visit change per donor and at equal donor weights | DONE | `test-visit-change.R` |
| 4 | Per-donor score draws exported (`cj_score_draws()`) | DONE | `test-scores.R`, `helper-scores-draws.R` |
| 5 | COMPASS adapter: FS, PFS, `Response()`, filters, heatmap | DONE | `test-compass.R` |
| 5 | COMPASS without a control sample: `missing_control = "donor"` or `"population"`, `control_source` | DONE | `test-compass.R`; calibration `scripts/75_calibrate-compass-missing-control.R`, `runs/2026-11-07_compass-missing-control/`; `docs/adr/0003` |
| 5 | Label-free calibration of non-paired COMPASS calls: `missing_control_calibrate = TRUE` | PARTIAL | `test-compass-calibrate.R`; `scripts/76` to `79`, `runs/2026-11-08_u2-*`, `runs/2026-11-09_u2r2-*`; `docs/adr/0003`. Clean held-out (seeds 41-50, 40 donors): FPR within 0.05 of paired in 12 of 12 cells. TPR bar 10 of 12 (two rare, weak cells at 25% missing). `"donor"` visit shift partly repaired by the switch to the same-visit population control. Needs 15 or more paired donors per run; costs 12 to 14 times the plain call. |
| 5 | COMPASS convergence diagnostics: `n_runs`, `diagnostics`, R-hat/ESS of gamma, FS, PFS, `unstable_call` | DONE | `test-compass-diagnostics.R`; validation `scripts/97_compass-convergence-sim.R`, `runs/2026-11-08_compass-diagnostics/` |
| 5 | COMPASS posterior differences | PARTIAL | not wrapped; `COMPASS::PosteriorDiff()` runs on `cj_compass()$fits` and is not tested here |
| 5 | `cj_mimosa2()`: native implementation of the MIMOSA2 design (ADR 0004), baseline-adjusted calls for one subset (`"any"` or one marker), four samples per donor; `missing_control = "marginal"` for one missing control | DONE | `test-mimosa2.R` (component log-likelihoods against hand sums, ordered means, q-values, counts by name, responders outrank non-responders on hand-built data, exclusions with reasons, marginal missing-control rows); native vs upstream in `runs/2026-11-04_mimosa2-parity/`; `"marginal"` calibrated in `runs/2026-11-04_mimosa2-marginal-{dev,heldout}/` (FPR bar 24 of 24 cells on held-out seeds 11-30, TPR loss 0.06 to 0.32) |
| 5 | `cj_baseline_report()`: COMPASS at baseline and follow-up beside the MIMOSA2 call, `agreement` class | DONE | `test-baseline-report.R` (stand-in COMPASS table reaching every class by name; not-assessed rows; one real COMPASS run) |
| 6 | Static ggplot2 plots | DONE | forest, response heatmap, contribution bars in `test-scores.R` |
| 7 | Adoption API: `cj_gate()` (tail, quantile, manual), `cj_plot_gate()`, `cj_from_cells()` | DONE | `test-adoption.R`; tail gate misclassifies under 0.5% of simulated cells |
| 7 | Readers: `cj_read_fcs()` (pre-gated FCS, spillover from the file), `cj_cells_from_gatingset()` | DONE (GatingSet part local only) | `test-read-cells.R` against flowCore by hand on flowWorkspaceData `CytoTrol` files and `gs_manual`; CI installs neither flowCore nor flowWorkspaceData, so the tests skip there |
| 7 | `cj_gate()` hardening: `method = "auto"`, `by` batch cuts, QC flags (`cut_in_broad_tail`, `high_control_positive`, `batch_cut_spread`), `cj_cofactor()` / `cofactor = "auto"`, `cj_simulate_gate_cells()` | DONE (`method = "auto"` is experimental: the held-out FlowJo comparison did not support it) | `test-gate-hardening.R`; dev evidence in `runs/2026-11-08_gating-hardening/dev*`; held-out FlowJo comparison in `heldout_*` |
| 7 | `cj_gate(method = "mindensity")` via `openCyto::gate_mindensity()` | DONE (local only) | `test-adoption.R` with `CJ_TEST_OPENCYTO=true`, openCyto 2.22.0 under R 4.5.1. The earlier crash was a host library clash (a user-library `XML` built against the system libxml2), not this package; CI does not install openCyto. `method = "tail"` is in-package: current openCyto has no `gate_tail()` |
| 7 | `cj_fluorescence_scores()`: NICR, Poly-NICR, frequency/intensity split, FIS, intensity-weighted FS/PFS (experimental), bootstrap | DONE | `test-adoption.R` (per-cell brute force, exact split, bootstrap); descriptive only, no Stan model, bootstrap covers within-sample cell sampling only |
| 7 | `cj_run_compass()` | DONE | `test-adoption.R` end to end with COMPASS 1.48 under R 4.5.1 |
| 7 | `cj_compass_report()`, `cj_plot_compass_report()`: COMPASS call beside C0-EB (default, empirical-Bayes pooled, no Stan), C0 or m1dm RMS/CPS with 90% intervals, NICR from cells, concordance flag | DONE | `test-compass-report.R`, `test-conjugate-eb.R` (COMPASS-free via a stand-in table, plus a real COMPASS run); simulation and RV144-dev validation in `runs/2026-11-08_compass-magnitudes/` |
| 7 | `cj_run_polyicsflow()` | DONE (local only) | `test-adoption.R` end to end with polyICSFlow 0.99.4 (SoegaardLab SHA 9dc2b37) under R 4.5.1; CI does not install polyICSFlow, so the test skips there |
| 7 | `cj_quickstart()` | DONE | `test-adoption.R` (without `fit`) |
| 6 | ComplexHeatmap output | PARTIAL | plots are ggplot2; ComplexHeatmap is only a Suggests entry and is not used |
| - | Pilot benchmark | PARTIAL | 13 of 21 scenarios x 30 seeds run; see "Pilot results" below |

## Known blockers

* `cj_from_gatingset()` has never run. flowWorkspace segfaults on the
  development host, even after its dependency chain was rebuilt from source.
* The polyICSFlow adapter is tested only on a hand-built table, and
  `cj_run_polyicsflow()` has not run through the installed package.
* openCyto cannot be loaded on the development host, so the `mindensity` gate is untested.
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

## Held-out generator D (2026-11-02): `m1dm` not vetoed

D was pre-registered before any candidate ran (outer repo commit `fed35d9`)
and opened only after the default was chosen. It has two parts:
- D1: a cell-level phenotype mixture;
- D2: generator A with t3 donor effects.

Both ran as rare_weak and moderate cells at J = 3 and 4, 30 seeds each (240
replicates, seed offset 900000). Source:
`runs/2026-11-02_redesign-D2-summary/calibration.tsv` and `false_sign.tsv`.
A first attempt (`runs/2026-11-02_redesign-D`) failed in the harness before
simulating; no model was fitted.

- **RMS:** `m1dm` covers 0.84 to 0.95 and passes all 8 cells within 2 SE.
  `m1` covers 0.78 to 0.92 and C0 covers 0.91 to 0.96.
- **CPS:**
  - `m1dm` covers 0.93 to 0.99. Three cells over-cover beyond the band
    (0.985 to 0.993); none under-cover.
  - `m1` under-covers on the phenotype mixture (D1: 0.41 to 0.48).
  - C0 fails D2 moderate (0.01 to 0.23).
- **CCS (experimental):** `m1dm` covers 0.94 to 0.997, all at or above
  nominal.
- **False signs:** `m1dm` 0, `m1` at most 0.002, C0 up to 0.99.
- **RMS RMSE:** `m1dm` and `m1` are 25 to 45% below C0 in D1 and in the
  D2 rare_weak cells, about 0.0004 vs 0.0006. In the D2 moderate cells all
  three are about equal (0.0011 to 0.0012).

Verdict: D does not veto `m1dm`. Its failures on D are conservative, the same
pattern as in A and C.

## Main limitation: CCS is not calibrated

**The generator-B RMS and CPS collapse is fixed by `m1dm`; CCS is still not
calibrated. CPS intervals are conservative when donors are homogeneous (coverage 0.98 to 1.00) and slightly narrow under marker-specific donor heterogeneity (0.81 to 0.84 of a nominal 0.90).** Under generator B the
pilot's M1 gave RMS coverage 0.23 and CPS 0.07. In the redesign benchmark,
`m1dm` gives RMS 0.90 and CPS 0.83 there (`calibration.tsv`, table above).

- **CCS** covers only 0.34 to 0.37 under `m1dm` and 0.27 to 0.39 under `m1`
  in generator B. No model covered it, so CCS is marked experimental: the
  `experimental` column, one message per `cj_scores()` call, and the
  "(experimental)" plot label. The C0 baseline gives CCS false signs up to
  0.995, so dropping the model is no fix.
- **CCS support gate:** CCS is withheld (`not_estimable`) below 20
  polyfunctional cells in the control or the stimulated sample, or above a
  posterior CV of 0.5 for the control's polyfunctional mass. In the 10-seed
  calibration (`runs/2026-11-03_ccs-gate/gate_chosen.tsv`, 24 donors, J = 3)
  no contrast passes in the three sparse regimes, where ungated generator-B
  coverage is 0.33. In dense regimes (5% any-positive) the gate passes 25% to
  100% of contrasts, with coverage 0.87 to 0.92 and false signs 0 to 0.013
  under generator B, and 0.99 to 1.00 (over-coverage) under generator A. Coverage
  is against the realised sample probabilities; against the structured ones it
  is 0.27 to 0.48 under B. The count condition decides; the CV condition adds
  nothing in these regimes and is a backstop. Regimes between the sparse and the
  dense ones were not simulated.
- **CPS under `m1dm`** covers 0.98 to 1.00 in generators A and C for a nominal
  0.90, so its intervals are wider than needed. In generator B it covers
  0.81 to 0.84. The cause of the over-coverage is unverified; the Dirichlet
  residual may inflate ratio-type scores.
- Generators A, B and C are simulations. Real data may differ from all three.

## Real-data results and limitations

COMPASS outperforms the default `m1dm` for responder classification on RV144 and TB.
The numbers come from the outer repo, [cytoJointR-reproducibility](https://github.com/mdmanurung/cytoJointR-reproducibility). Details and sources are in `runs/2026-11-03_realdata/README.md`; the benchmark is in `benchmark/README.md`.

| Data set | AUC, `m1dm` P(RMS > eps), sample level | AUC, COMPASS FS | Donors |
|---|---|---|---|
| RV144, vaccine against placebo | 0.755 | 0.977 | 226 vaccine, 36 placebo |
| TB_s, QFT+ against QFT- | 0.758 | 0.987 | 22 and 18 |
| TB_u, QFT+ against QFT- | 0.649 | 0.864 | 22 and 18 |

- **False positives.** The sample-level rule P(RMS > 1e-4) > 0.95 calls 8 of 36 RV144 placebos positive. COMPASS calls 0 of 36.
- **Structured level.** It ranks better on RV144 (AUC 0.932) and TB_s (0.917). It does not on TB_u (0.610). At the 0.95 rule it calls 2 of 22 QFT+ donors on TB_s and 0 on TB_u (TPR 0.091 and 0.000). At eps 1e-4 it still calls 5 of 36 RV144 placebos, and it shows no infection signal (odds ratio per SD 0.997).
- **Mean misfit.** The fitted structured mean exceeds the observed rates. For cells with 3 or more positive markers the ratio is 4.5 to 24 in RV144. The cause is not established.
- **No usable false-positive rate elsewhere.** HVTN078 has 5 placebos, so its AUC is uninformative. The SARS-CoV-2 donors had prior COVID-19, so Spike responses at baseline are expected.
- **SimpleCOMPASS** (21 donors): the intervals overlap and no method is clearly better.
- **Experimental responder-mixture model.** Branch `feat/donor-effects` (not merged) reaches RV144-dev AUC 0.98, against 0.974 for COMPASS. The paired difference is not distinguishable from 0 (best case +0.007, 95% CI -0.007 to 0.025). It stays at TB_u AUC 0.61 to 0.64, against 0.864 for COMPASS, and fails some simulation coverage guards. It does not beat COMPASS on the frozen criteria.

## Parameterisation

Default is noncentered. `centered = TRUE` is a speed option: it funnels when a
donor variance is near zero (pilot: R-hat up to 1.4 for null truths; 1.23 when
`donor_visit = TRUE` meets data without donor-by-visit variation). In the pilot
noncentered fits at 100 donors (median 503-643 s) were not slower than centered
ones (1020-1322 s, different scenarios).

## Deferred (from SCOPE)

A fluorescence model (M2: truncated-Gaussian conditional intensities with
posterior uncertainty for NICR, Poly-NICR, FIS; descriptive versions are
implemented, see item 7); RTMB engine; group terms
and `cj_group_contrast()`; donor-specific induced change across visits;
per-visit interaction change; continuous-time visits; Shiny app; hierarchical
missing-control prediction; CRAN/Bioconductor build layer.
