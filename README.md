# cytoJointR

Donor-paired Bayesian inference on polyfunctional intracellular cytokine
staining (ICS) data. The package models exact marker patterns within donors,
relative to each donor's matched control sample, and reports signed net scores
with posterior probabilities. Terms are defined in [GLOSSARY.md](GLOSSARY.md).

## Install

```r
# install.packages("remotes")
remotes::install_github("mdmanurung/cytoJointR", build_vignettes = TRUE)
# or, once published: install.packages("cytoJointR", repos = c("https://mdmanurung.r-universe.dev", getOption("repos")))

# Stan models compile on first use
cmdstanr::install_cmdstan()
```

`cmdstanr` comes from `https://stan-dev.r-universe.dev`. COMPASS (Bioconductor)
is optional and needed only for `cj_compass()` and `cj_run_compass()`; openCyto
(`cj_gate(method = "mindensity")`) and polyICSFlow (`cj_run_polyicsflow()`) are optional too.

## Example

```r
library(cytoJointR)
sim <- cj_simulate_counts(n_donors = 12, markers = c("IFNG", "TNF", "IL2"),
                          visits = c("v1", "v2"), N = 1500, donor_visit = TRUE, seed = 2)
summary(sim$data)
cj_plot_degree_bars(sim$data, "at_least")
fit <- cj_fit(sim$data, centered = TRUE, adapt_delta = 0.95, seed = 1)
eps <- c(RMS = 0.005, poly_absolute = 0.005, poly_marginal = 0.005, CPS = 0.02,
         CCS = 0.005, breadth = 0.05, marker_magnitude = 0.005)
sc <- cj_scores(fit, eps)
cj_plot_forest(sc, score = c("RMS", "CCS"))
```

From a cell data frame (one row per event) to scores:

```r
samples <- data.frame(sample_id = c("D1_ctrl", "D1_cmv"), donor_id = "D1", visit = factor("v1"),
                      parent_population = "CD4", condition = c("unstim", "CMV"))
gate <- cj_gate(cells, c("IFNG", "TNF", "IL2"), samples, transform = "asinh", cofactor = 150,
                control = "unstim")             # one pooled threshold per marker
cj_plot_gate(gate)                              # check every panel
x  <- cj_from_cells(cells, samples = samples, gate = gate)
fs <- cj_fluorescence_scores(x, n_boot = 200, seed = 1)  # NICR, Poly-NICR, FIS, iw_FS/PFS
cp <- cj_run_compass(x)                         # needs COMPASS
rp <- cj_compass_report(x)                      # COMPASS call + RMS with 90% interval, concordance
cj_plot_compass_report(rp)
```

`cj_quickstart()` runs these in one call. `vignette("from-cells")` has the details.

`vignette("cytoJointR")` walks through the full workflow.
`vignette("comparison")` maps COMPASS and polyICSFlow features to this package.

## Scope status

**Limitation:** the default model, `m1dm`, covers RMS well, including when donors
vary in marker-specific ways (2026-11-02 benchmark). Its CPS intervals are
conservative when donors are homogeneous and slightly narrow (0.81 to 0.84) under
marker-specific donor heterogeneity. CCS is experimental: no model covered it
under donor heterogeneity (coverage 0.27 to 0.39 of a nominal 0.90 in generator B), so do not
treat CCS intervals or signs as calibrated. The model default is a post hoc
choice; see `docs/adr/0002-default-model-m1dm.md`.


This is a development version (v0.1 in progress). The input layer,
descriptive layer, count model (`m1dm`, with `m1` for ablation), score layer, visit change and the COMPASS
adapter are implemented and tested. `cj_from_gatingset()` has never run, the
pilot benchmark is not finished, and CI has never run. See
[IMPLEMENTATION_STATUS.md](IMPLEMENTATION_STATUS.md) for each item and its
evidence, and [docs/SCOPE.md](docs/SCOPE.md) for what is in and out of v0.1.

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
