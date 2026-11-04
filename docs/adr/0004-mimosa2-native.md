# 0004. A native MIMOSA2 model, written from its design

Date: 2026-11-04. Status: accepted.

## Context

MIMOSA2 (github.com/RGLab/MIMOSA2, version 0.99.0.900, GPL-3) gives
baseline-adjusted positivity calls for intracellular cytokine staining. It has
no paper; its design is known from its README and source. cytoJointR is MIT.
Copying MIMOSA2's code would make cytoJointR GPL-3. The user chose to keep MIT
and implement the model afresh from its design (2026-11-04).

Two needs the upstream package cannot meet: donors with a missing control
sample (its initialisation and one-sided rule divide by every sample's total),
and an install from a mainstream repository.

**Licence condition.** The rewrite reimplements the design of an unpublished
method, and the only description of that design is GPL-3 code. Implementing it
from the description below, in new code, is our basis for MIT. For certainty,
the alternatives are GPL-3 for cytoJointR or the author's permission.

## The model

One fit covers the donors of one stimulation, parent population and follow-up
visit. Each donor has four counts of positive cells out of a total: stimulated
at follow-up (`s1`), control at follow-up (`u1`), stimulated at baseline
(`s0`) and control at baseline (`u0`).

**Levels.** Four proportion levels `k = 1..4`, each a beta distribution with
mean `mu_k` and precision `exp(nu_k)`, ordered `mu_1 > mu_2 > mu_3 > mu_4`.
A donor's samples that share a level in a component share one proportion drawn
from that beta: their counts and totals are summed into one beta-binomial term.

**Components.** Eleven components say which samples share which level.
Responder components are 1 to 4.

| Component | Label | Terms (samples pooled in braces) |
|---|---|---|
| 1 | R1 | `s1` at 1, `s0` at 2, `u1` at 3, `u0` at 4 |
| 2 | R2 | `s1` at 1, `u1` at 3, `{s0, u0}` at 4 |
| 3 | R3 | `{s1, s0}` at 1, `u0` at 3, `u1` at 4 |
| 4 | R4 | `s1` at 1, `s0` at 2, `{u1, u0}` at 4 |
| 5 | NR1 | `{s1, u1}` at 3, `{s0, u0}` at 4 |
| 6 | NR2 | `{s1, u1}` at 3, `s0` at 2, `u0` at 4 |
| 7 | NR3 | `{s1, u1, s0, u0}` at 4 |
| 8 | NSR1 | `{u1, u0}` at 1, `{s1, s0}` at 4 |
| 9 | NSR2 | `u1` at 1, `{s1, s0, u0}` at 4 |
| 10 | NSR3 | `u0` at 1, `{s1, u1, s0}` at 4 |
| 11 | NSR4 | `{s1, s0}` at 2, `{u1, u0}` at 4 |

**One-sided rule.** A donor whose observed
`delta = s1/N_s1 - u1/N_u1 - s0/N_s0 + u0/N_u0` is negative cannot be a
responder: its responder components get likelihood 0.

**Priors (MIMOSA2's default values).** `nu_k ~ Gamma(shape 11/4, rate 1/2)` on
each log precision. Normal penalties with mean `1e-4` and SD `5e-4` on five
mean contrasts: `mu_1 - mu_3 - mu_2 + mu_4`, `mu_1 - mu_3`, `mu_2 - mu_4`,
`mu_1 - mu_2` and `mu_3 - mu_4`.

**Fit.** EM to the posterior mode. The E-step gives each donor's component
probabilities `z` from `pi_k` times the component likelihood. The M-step sets
`pi` to the column means of `z` and maximises the expected complete-data log
posterior over the eight level parameters.

**Output.** `P(response)` is the sum of `z` over components 1 to 4. The q-value
is that of Newton et al. (2004): order donors by `P(response)`, decreasing, and
take the running mean of `1 - P(response)`.

## Our choices (where this differs from upstream MIMOSA2)

- **Mixing weights** enter as `log(pi_k)`. Upstream uses `log1p(pi_k)`.
- **Each level has its own precision.** Upstream's R1 gives `s0` the mean of
  level 2 with the precision of level 1.
- **Ordering by construction.** `logit(mu_4) = a`,
  `logit(mu_k) = logit(mu_{k+1}) + exp(b_k)`. Upstream sorts the means after
  each step, which can swap a mean away from its precision.
- **The one-sided rule is hard** (likelihood 0). Upstream sets those entries to
  1 below the smallest log-likelihood of the initial fit.
- **Optimiser:** `stats::optim(method = "BFGS")` on the unconstrained
  parameters. Upstream uses optimx with newuoa and bobyqa.
- **Initialisation** from the data: donors with `delta > 0` and `s1/N_s1 >
  u1/N_u1` start in the responder components (spread evenly), the others in
  the non-responder and non-specific components. Level means start at the
  responders' mean `s1` and `s0` proportions and everyone's mean `u1` and `u0`
  proportions, then sorted. Log precisions start at the prior mode, 3.5.
- **Convergence:** stop when the objective changes by less than `tol`
  (relative), or after `maxit` iterations. The fit reports whether it
  converged.

## Missing samples (stage 2, planned before any result)

A missing sample has total 0 and count 0. Its beta-binomial term is then
exactly 1, so it drops out of the likelihood: this is marginalising the missing
sample, valid when it is missing at random.

- A missing **stimulated** sample removes the information on response: such
  donors stay excluded.
- A donor with **both controls missing** stays excluded.
- **One control missing:** the donor enters. The initialisation and the
  one-sided rule use the donor's other control in place of the missing one.
  This assumes the background is stable across visits, the same assumption as
  `cj_compass(missing_control = "donor")`.
- **What it costs:** with `u1` missing, telling responders from non-responders
  rests on whether `s1` sits at level 1 or level 3, an absolute level rather
  than a within-donor contrast. The calibration simulation therefore includes
  between-donor background heterogeneity and a background shift between
  visits. It uses dev seeds for development, then one clean test on held-out
  seeds, as for `missing_control` in `cj_compass()`.

## Validation of the rewrite

1. The component log-likelihoods match upstream `cll()` on random parameters
   to 1e-10 for components 2 to 11. Component 1 matches when levels 1 and 2
   share a precision (our precision fix).
2. On upstream's `simulate_MIMOSA2()` data and the package's hand-built data,
   the native fit and upstream agree on calls, with Spearman correlation of
   `P(response)` reported. The numbers are in a run README of the
   reproducibility repository.
