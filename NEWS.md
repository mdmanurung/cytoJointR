# cytoJointR 0.0.0.9000

* Package started from the unified design bundle of 2026-11-01. Reference
  utilities (`cj_bits`, `cj_compress`, likelihood, IPF counterfactual,
  `cj_score_reference`, `cj_summarize_signed`) carried over unchanged; they
  passed 35/35 inherited expectations under R 4.5.1.
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
