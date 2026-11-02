# cytoJointR 0.0.0.9000

* Package started from the unified design bundle of 2026-11-01. Reference
  utilities (`cj_bits`, `cj_compress`, likelihood, IPF counterfactual,
  `cj_score_reference`, `cj_summarize_signed`) carried over unchanged; they
  passed 35/35 inherited expectations under R 4.5.1.
* Model redesign (2026-11-02): `cj_fit()` defaults to `model = "m1dm"`, the M1
  mean with a Dirichlet-multinomial sample residual; `model = "m1"` stays for
  ablation. The default is a post hoc choice after the pre-registered rule
  found no calibrated model. CPS intervals are conservative when donors are homogeneous (coverage 0.98 to 1.00) and slightly narrow under marker-specific donor heterogeneity (0.81 to 0.84 of a nominal 0.90). See
  `docs/adr/0002-default-model-m1dm.md`.
* CCS is experimental: `cj_scores()`, `cj_score_draws()` and `cj_patterns()` add
  a logical `experimental` column (`TRUE` for CCS), `cj_scores()` emits one
  message per call, and plots label it "CCS (experimental)".
* CCS support gate: `cj_scores()`, `cj_score_draws()` and `cj_patterns()` gain
  `ccs_gate` (default `list(min_poly_cells = 20, max_cv = 0.5)`; `NULL` switches
  it off). A contrast keeps its CCS rows but gets `status = "not_estimable"`, a
  `reason` and `NA` values unless the control and the stimulated sample each hold
  at least `min_poly_cells` polyfunctional cells and the posterior CV of the
  control's polyfunctional mass is at most `max_cv`. A `cj_draws` has no counts,
  so its CCS rows get `status = "support_unknown"`; pass a `cj_fit` or
  `ccs_gate = NULL`. A `reason` column is added to all three outputs. CCS stays
  `experimental`. `cj_visit_change()` is not gated.
* `cj_scores_from_probs(min_functional=)` is deprecated and ignored.
* RTMB removed from Suggests: Stan via cmdstanr is the only engine in v0.1.
* Input layer: `cj_from_counts()`, `cj_from_events()`, `cj_from_gatingset()`
  (not yet exercised: flowWorkspace fails to load on the development host),
  `cj_from_polyicsflow()`. Contrasts are derived per donor x visit x parent
  population; a missing control sample makes that contrast's net scores
  unavailable.
* Descriptive layer: `cj_describe()` with signed net percentages and a
  display-only clipped column; pattern heatmap, degree bars, pattern scatters.
* Score layer: `cj_score_draws()`, `cj_scores()`, `cj_patterns()`,
  `cj_markers()`, `cj_visit_change()`; forest plot, response heatmap,
  contribution bars. IPF failures are flagged per draw, never dropped.
* `cj_compass()`: runs native COMPASS per stimulation x visit x parent
  population and returns the native `COMPASSResult` objects.
* `cj_simulate_counts()` now defaults to `donor_visit = TRUE`, matching
  `cj_fit()`; the old mismatch gave R-hat 1.09 when users copied the pattern.
* Pilot benchmark (2026-11-02): good recovery when data come from M1, but
  coverage collapses under logistic-normal donor heterogeneity (RMS 0.23).
  See IMPLEMENTATION_STATUS.md, "Main limitation".
* `cj_design()`, `cj_fit()` and `cj_prior_predict()` stop when the samples hold
  more than one parent population. New `cj_subset()` keeps one.
* `cj_fit()` writes each fit's CmdStan CSVs to its own directory. Before, two
  fits with the same `set.seed()` in the same minute got identical file names,
  and garbage collection of the first fit deleted the second fit's draws.
* `cj_fit(model = "m1dm")` (new default) adds a Dirichlet-multinomial sample
  residual to M1: `y_s ~ DM(phi * softmax(eta_s))`, `inv_phi = 1 / phi ~
  half-normal(0, 0.05)`. `model = "m1"` is the old model; its draws are
  identical at a fixed seed. `level = "structured"` returns `softmax(eta_s)` in
  `fit$draws`; `fit$draws_sample` and `fit$draws_structured` hold both. The
  default model may change after the planned model comparison.
* `cj_simulate_counts(truth = list(phi = ))` simulates the same residual (default
  `Inf`, the old output) and returns `mu` (structured) beside `p` (realised).
  `cj_prior_predict(model = )` draws `inv_phi` for `"m1dm"`. `cj_design()` gains
  `overdispersed`; `cj_default_priors()` gains `inv_phi_sd`.
* The compiled Stan executable is cached under a name that includes a hash of
  `m1.stan`, so checkouts with different Stan code no longer overwrite each other.
* `cj_compass()` now computes the COMPASS posterior, so `COMPASS::PosteriorDiff()`
  and its relatives work on `$fits`.
* `cj_describe()` no longer errors on data without contrasts.
* Internal `cj_fit_dirichlet()`: independent conjugate `Dirichlet(a0 + y_s)` per
  sample, no Stan; the baseline the hierarchical models must beat.
