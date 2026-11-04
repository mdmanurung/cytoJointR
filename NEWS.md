# cytoJointR 0.0.0.9000

* `cj_read_fcs()` reads pre-gated FCS files (one per sample) into the `cells` table, applying
  each file's spillover matrix (`compensate = TRUE`); markers are chosen by channel name or
  exact channel description. `cj_cells_from_gatingset()` does the same for one population of a
  flowWorkspace `GatingSet`. Both feed `cj_gate()`, `cj_from_cells()` and `cj_quickstart()`.

* `cj_baseline_report()`: one row per donor and follow-up visit with COMPASS at the baseline and
  follow-up visits beside the baseline-adjusted `cj_mimosa2()` call, and an `agreement` class
  (`both`, `pre_existing`, `compass_only`, `mimosa2_only`, `neither`, `not_assessed`).
  `pre_existing` marks donors COMPASS calls at both visits whose response MIMOSA2 does not
  see rise.

* `cj_mimosa2()` is a native implementation of the MIMOSA2 design (ADR 0004; no upstream code, no
  extra package) for baseline-adjusted positivity calls: one fit per stimulation, parent
  population and follow-up visit, with each donor's stimulated and control samples at the baseline
  and the follow-up visit. One subset per fit: `positive = "any"` (any marker positive) or one
  marker. Donors without all four samples are listed in `excluded`. Scores give the four counts,
  `delta`, `mimosa2_p_response`, `mimosa2_fdr` (q-values within one fit), `mimosa2_response` and
  `control_source`; `cj_mimosa2_fit` objects are in `fits`. `missing_control = "marginal"` lets a
  donor with one missing control enter. In a held-out simulation its false-positive rate stayed
  within 0.05 of paired donors in 24 of 24 cells; it costs power (TPR 0.06 to 0.32 lower).

* `cj_compass(n_runs = , diagnostics = )` checks convergence of the COMPASS call. COMPASS's
  `replications` are passes of one chain, not independent chains, so the check refits with seeds
  `seed`, `seed + 1`, ... (`diagnostics = TRUE` means at least 5 runs, so a call never ties) and treats the fits as chains.
  Per donor and stimulation, `scores` gains rank-normalised R-hat and bulk ESS of the posterior
  response probability (mean of `gamma`), FS and PFS (`compass_gamma_rhat`, `compass_FS_ess`, ...), the
  number of runs that call the donor (max mean `gamma` > 0.5), `unstable_call` (the runs disagree on
  that call; with `rhat_max = 1.05` also when an R-hat exceeds it) and `compass_call_unsupported` (at
  most half of the runs make the reported call). R-hat does not flag by default: it found no wrong
  call that run disagreement missed and flagged 27% of simulated donors.
  Simulation behind the defaults (40 donors, 4 markers, 12 weak-effect data sets, a 10x-iterations
  reference from 4 seeds; 18 of 480 donors had a wrong reported call at default iterations): with 5 runs
  `unstable_call` has recall 1.00 and precision 0.35, `compass_call_unsupported` recall 1.00 and
  precision 0.95; with 4 runs the latter's precision is 0.50 (ties). At 0.1x iterations their recalls are
  0.95 and 0.63. R-hat above 1.05 flagged 27% of donors (precision 0.12). On 131 RV144 dev donors with
  5 runs, 1 donor was unstable at default and at 2x iterations, none unsupported. `cj_compass()$diagnostics$summary` counts flagged donors and
  gives the runtime. The reported scores are those of the first run, so `n_runs = 1` (default) is
  unchanged. `cj_run_compass()` passes both arguments on. Needs the posterior package. The check costs
  about `n_runs` times one fit.

* `cj_compass_report()`: one table per donor, stimulation and visit with the COMPASS
  call (FS, PFS, response probability, `control_source`) beside cytoJointR
  magnitudes with 90% intervals (RMS, CPS, per-marker magnitude) and, from cells,
  NICR and Poly-NICR with a cell bootstrap. The default magnitude is C0-EB,
  an empirical-Bayes pooled conjugate model that needs no Stan (`"c0"` is the
  plain independent Dirichlet); `magnitude_model = "m1dm"` uses the
  structured posterior. C0-EB warns for a parent population with fewer than 8 donor-visit blocks, fills a partial `epsilon` from the defaults, and warns when `a0`, `n_draws` or `fit_args` are ignored. A `concordance` column flags a COMPASS call whose RMS
  interval includes 0, and the reverse. `cj_plot_compass_report()` draws the
  table as a forest plot with the COMPASS response as point colour.

* Hardened `cj_gate()` (all additions; defaults unchanged).
  `by = "<samples column>"` finds one cut per batch on that batch's control events
  and `cj_from_cells()` applies it (events are moved by the batch cut minus the
  pooled cut); batches with too few control events (`min_batch_events`, default
  1000) warn and take the pooled cut. `method = "auto"` chooses tail or quantile
  per marker from a label-free shape statistic of the control events (the tail
  ratio, with a cap on the control mass above the tail cut), with constants fixed
  from simulation. `$qc` and `print()` carry three flags per marker:
  `cut_in_broad_tail`, `high_control_positive`, `batch_cut_spread`;
  `cj_plot_gate()` titles each panel with them and draws the batch cuts.
  `cj_cofactor()` and `cofactor = "auto"` give a label-free asinh cofactor from the
  width of the negative peak; a cofactor named by marker is accepted.
  `cj_simulate_gate_cells()` simulates controls with clean and heavy upper tails,
  weak positives in the control and batch offsets, with a truth column. The rule,
  its dev evidence and the held-out FlowJo comparison are in the vignette
  `from-cells` and `runs/2026-11-08_gating-hardening/`. On two public FCS datasets
  (9 markers), the old tail gate had mean kappa 0.604 against the expert gates;
  tail + `by` + `cofactor = "auto"` had 0.825. `method = "auto"` did not help there
  (0.761) and is experimental.

* `cj_compass(missing_control = c("drop", "donor", "population"))` runs COMPASS on
  stimulated samples that have no control sample. `"drop"` (default) is the old
  behaviour. `"donor"` borrows the same donor's control from the nearest other
  visit (ties to the earlier visit; no control at all falls back to
  `"population"`, flagged). `"population"` passes an all-zero control row, which
  COMPASS compares with the other donors' control distribution. Scores gain
  `control_source` and `control_fallback`, a message counts the unpaired rows,
  and `cj_run_compass()` passes the option on. The assumptions and the
  simulation evidence are in `docs/adr/0003-compass-missing-control.md`.

* Adoption API: from a cell data frame to scores in a few calls.
  `cj_gate()` finds one threshold per marker on all control events pooled
  (`"tail"` default, `"quantile"`, `"mindensity"` through openCyto, or `"manual"`)
  on an explicit `transform` (never assumed); `cj_plot_gate()` shows the pooled
  control density with the threshold. `cj_from_cells()` builds the `cj_data` and
  keeps per-pattern sums of `X - t` and `(X - t)^2`. `cj_run_compass()` and
  `cj_run_polyicsflow()` return tidy per-donor tables. `cj_quickstart()` chains
  them. `method = "tail"` is implemented in the package: current openCyto has no
  `gate_tail()`. openCyto, COMPASS and polyICSFlow are `Suggests`.
* `cj_fluorescence_scores()`: descriptive, model-free NICR, Poly-NICR, FIS and
  an exact frequency/intensity split per donor, stimulation, visit and marker,
  with a cell-level bootstrap. Experimental `iw_fs` and `iw_pfs` apply COMPASS's
  FS and PFS weights to the net burden instead of the response probability.
  Definitions: `?cj_fluorescence_scores` and `GLOSSARY.md`. `cj_subset()` keeps
  the intensity summaries.

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
