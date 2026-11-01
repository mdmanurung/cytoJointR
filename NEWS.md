# cytoJointR (development version)

* Package started from the unified design bundle of 2026-11-01. Reference
  utilities (`cj_bits`, `cj_compress`, likelihood, IPF counterfactual,
  `cj_score_reference`, `cj_summarize_signed`) carried over unchanged; they
  passed 35/35 inherited expectations under R 4.5.1.
* `cj_scores_from_probs(min_functional=)` is deprecated and ignored.
* RTMB removed from Suggests: Stan via cmdstanr is the only engine in v0.1.
