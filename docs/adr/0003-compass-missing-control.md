# 0003. COMPASS without a control sample runs in the adapter

Date: 2026-11-03, calibration added 2026-11-08. Status: accepted.

## Context

`COMPASS::COMPASS()` and `COMPASSContainer` need a paired unstimulated sample
for every stimulated sample. Real studies lose controls: a failed aliquot, a
shared control that never arrived. `cj_compass()` dropped those contrasts.

COMPASS integrates out each donor's control proportions `p_u ~ Dir(alpha_u)`.
A control row of zero counts therefore adds nothing to the likelihood, and that
donor's stimulated counts are compared with the population control
distribution that the other donors define. The exported
`COMPASS::SimpleCOMPASS()` accepts such a row. The pairing is enforced only by
COMPASS's input layer and by our adapter.

## Decision

`cj_compass(missing_control = c("drop", "donor", "population"))`:

- `"drop"` (default): unchanged, the contrast is excluded.
- `"donor"`: the same donor's control from the nearest other visit (same parent
  population; ties to the earlier visit). A donor with no control falls back to
  `"population"`, flagged.
- `"population"`: an all-zero control row.

Every output row says where its control came from (`control_source`). The
change lives in the adapter. We call only exported COMPASS functions
(`SimpleCOMPASS()`, `FunctionalityScore()`, `Response()`, `PosteriorDiff()`) and
do not fork or patch COMPASS (Artistic-2.0, and a fork would carry
maintenance and a second copy of the sampler). `compute_posterior` is the same
unexported helper the adapter already used.

## Assumption

`"population"` assumes the donor's unstimulated pattern distribution looks like
the other donors' controls. A donor whose background differs (high background,
a different assay day) gets a wrong comparison: the excess background reads as
a response. `"donor"` assumes the background does not change between visits. A
real visit effect on the background biases the borrowed contrast. The size of
both effects on simulated data is in `runs/2026-11-07_compass-missing-control/`
of the reproducibility repo. Neither option is the default.

## Caveats built in

Both options are experimental: a message says so and gives the simulated FPR.
`"donor"` compares the donor's available controls across visits (any-positive
fraction; proportion test p < 0.01 and a 1.25-fold difference, fixed a priori)
and sets `control_shift_check` to `ok`, `differs` or `untestable`, warning on the
last two.

## Consequences

- Default behaviour and outputs are unchanged apart from two new columns.
- A run in which no donor has a control cannot borrow a population control; it
  is excluded with a reason.
- COMPASS quantities for a zero-control row (FS, PFS, response probabilities,
  `PosteriorDiff()`) are finite and are relative to the population control.
- If a COMPASS release stops accepting a zero row, `"population"` fails with
  COMPASS's message and the group is listed as failed in `groups`.

Rejected: forking COMPASS; imputing a control from the population mean counts
(it invents data and gives the donor control a fake precision).

## Calibration of non-paired calls (2026-11-08)

`cj_compass(missing_control_calibrate = TRUE)` addresses the liberal call: the
`gamma > 0.5` rule gave FPR 0.19 to 0.39 (`"population"`) among non-paired rows against
0.05 to 0.10 for paired rows.

Method (label-free). In each COMPASS run, paired donors are pseudo-removed in folds
of 5, up to `max_pseudo` (default 40), and the run is refitted with the same
`missing_control` rule. A pseudo donor has a gamma with its control missing and a gamma
with its own control. Isotonic regression of the paired gamma on the pseudo-missing
gamma maps the non-paired rows, per run and per kind of control (borrowed or
population). The calibrated gamma never exceeds the raw gamma. A pseudo fit has its pseudo-removed donors
plus the truly missing ones with no control, so the pseudo fit sees a higher missing fraction
than the real fit (by up to 5 donors per fold); its gamma is slightly more liberal than the
real missing-control gamma would be, which biases the map toward being too lenient at high
missing rates. The warning on a poor pseudo AUC applies only to runs that used the map. The call is
`compass_call` (calibrated gamma > 0.5). No responder labels are used.

Choices and why (simulation of the 2026-11-07 design, 4 markers, 50% responders;
dev seeds 1 to 10, held-out seeds 11 to 20; `runs/2026-11-08_u2-compass-calibrate*`):

- Isotonic beat the quantile map (quantile max FPR gap +0.20 to +0.21 on dev, isotonic
  +0.03 to +0.07). An empirical-null cut (95th percentile of the pseudo-missing gamma of
  paired non-calls) did no better than isotonic.
- A bigger pseudo set helps: `"population"` max FPR gap over paired +0.107 with 10
  pseudo donors, +0.065 with 20, +0.010 with 40. Default 40.
- Minimum: 15 paired donors per run (at least 3 pseudo donors above and 3 not above
  gamma 0.5). With 10 to 14 pseudo donors the map gave FPR +0.11 over paired. Below the
  minimum the run uses the fallback cut: raw gamma 0.85 sent to 0.5 (Option 2). The cut
  is the point where mean non-paired minus paired FPR crossed zero on dev seeds 1 to 10
  of the 2026-11-07 run. It depends on that simulation design. It is safe on FPR and
  costs power at a rare, weak effect.
- Tightening only (calibrated call never above the raw call) removed the extra false
  calls that the map added at a rare, weak effect (held-out FPR gap +0.056 became
  +0.036 at worst).

Held-out result (seeds 11 to 20, frozen candidate, 40 donors; non-paired against paired
rows of the same fit):

| Setting | FPR raw | FPR calibrated | FPR paired |
|---|---|---|---|
| `"population"`, moderate | 0.21 to 0.24 | 0.06 to 0.09 | 0.07 |
| `"donor"`, moderate, stable background | 0.05 to 0.14 | 0.03 to 0.05 | 0.07 |
| rare, weak (both) | 0.02 to 0.06 | 0.01 to 0.06 | 0.005 to 0.02 |

Card bars: non-paired FPR within paired + 0.05 held in 12 of 12 cells (largest gap
+0.037); non-paired AUC fell at most 0.008. Non-paired TPR fell by at most 0.05 in
8 of 12 cells. The failures are `"donor"` at the rare, weak effect (TPR -0.06, -0.08,
-0.19 at 10, 25, 50% missing: with half the controls missing, fewer than 15 paired
donors have a usable borrowed control, so the run uses the fallback cut) and
`"population"` moderate at 10% missing (-0.05, from 1.00 to 0.95).

Where it does not hold, in plain terms. With `"donor"` and a visit effect, calibration
does not repair the false-positive rate at the visit whose control is borrowed. With two
visits `control_shift_check` returns `"untestable"` for every borrowed row, so it cannot
tell a shifted from an unshifted background. The detector to watch is the pseudo-AUC
warning (`$calibration$table$pseudo_auc` below 0.8). Details: `"donor"` with a visit effect (background x1.5 at visit 2): the
later visit keeps FPR 0.33 to 0.95 after calibration (0.82 to 0.95 raw at 10 to 25%
missing at the moderate effect; 0.57 against 0.93 at 50%). Pseudo-missing gamma cannot
rank donors there, and the calibration says so: the pseudo AUC is below 0.8 in 100% of
those runs, and in 0% of the moderate-effect runs without a visit effect. At the rare, weak
effect the same warning also fires in 10 to 18% of runs without a shift, because few
paired donors have a call to compare. `control_shift_check` did not help in this design:
with two visits the donor has one other control, so it is always `"untestable"` (it
warns for every borrowed row, shifted or not) and never `"differs"`. So the card's
80% recall holds only in the sense that all rows warn, and the pseudo-AUC warning is the
usable detector. Controls missing not at random, pattern-specific visit effects and
more than one stimulation are not simulated.

Cost. 40 donors, 5000 iterations x 4 replications: the plain call takes 19 to 47 s and
the calibrated call 97 to 394 s (5 to 9 times; held-out median). `max_pseudo` caps
it: 20 costs about 5 times, 40 about 5 to 9 times.


Fallback ablation (dev seeds 1 to 10, before the clean held-out test). Keeping the raw gamma
for runs with too few pseudo donors, instead of the 0.85 cut, meets the FPR bar in 9 of 12
cells (largest gap +0.118, `"donor"` moderate) against 12 of 12 (+0.049) with the cut.
The cut costs power: TPR bar 10 of 12 against 12 of 12. The fallback stays because it earns
its place on the FPR bar; its constant remains simulation-derived.

Clean held-out test (seeds 21 to 30, scored once, package 1345d2d, same grid and 40 donors,
`runs/2026-11-08_u2-compass-calibrate-heldout-clean`). The earlier held-out seeds 11 to 20
were used twice and informed the tighten-only rule, so they are no longer a clean test; this is.
- FPR bar (non-paired <= paired + 0.05): 11 of 12 cells. The miss is `"population"`,
  moderate effect, 50% missing: 0.158 against 0.071 paired (+0.086). Seeds 26 and 28 keep
  almost all raw false calls there, with only about 20 paired donors for the map. Other cells
  have gaps of -0.017 to +0.017. The raw rule met the bar in 6 of 12 (largest gap +0.285).
- TPR bar (loss <= 0.05): 9 of 12. The three `"donor"` rare-weak cells fail (-0.07, -0.10,
  -0.24), through the fallback cut, as before.
- Non-paired AUC fell at most 0.004.
- Visit effect with `"donor"`: unchanged. The later visit keeps FPR 0.52 to 0.91 (moderate)
  and 0.55 to 0.72 (rare, weak). `control_shift_check` was `"untestable"` for every borrowed
  row, and the pseudo-AUC warning fired in 100% of those later-visit runs against 0% of moderate
  runs without a shift (14 to 23% at the rare, weak effect without a shift).
So the claim is: calibration brings non-paired FPR to within 0.05 of paired in most cells, not
all, at 40 donors with about 20 or more paired donors per run.
