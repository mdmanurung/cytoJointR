# 0003. COMPASS without a control sample runs in the adapter

Date: 2026-11-03. Status: accepted.

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
