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
is optional and needed only for `cj_compass()`.

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

`vignette("cytoJointR")` walks through the full workflow.
`vignette("comparison")` maps COMPASS and polyICSFlow features to this package.

## Scope status

**Limitation:** in the 2026-11-02 pilot, intervals were well calibrated when
data came from the model itself, but far too narrow (RMS coverage 0.23 of a
nominal 0.90) when donors varied in marker-specific ways. Treat intervals on
real data as optimistic until that is fixed.


This is a development version (v0.1 in progress). The input layer,
descriptive layer, M1 count model, score layer, visit change and the COMPASS
adapter are implemented and tested. `cj_from_gatingset()` has never run, the
pilot benchmark is not finished, and CI has never run. See
[IMPLEMENTATION_STATUS.md](IMPLEMENTATION_STATUS.md) for each item and its
evidence, and [docs/SCOPE.md](docs/SCOPE.md) for what is in and out of v0.1.
