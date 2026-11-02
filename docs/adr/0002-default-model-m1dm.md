# 0002. Default count model is m1dm, chosen post hoc

Date: 2026-11-02. Status: accepted. The decision is post hoc: the
pre-registered rule returned no default, and the user chose one anyway.

## Context

In the pilot, M1 covered well when data came from M1 itself but collapsed under
generator B (donor, stimulation and sample noise at pattern level): RMS
coverage 0.23, CPS 0.07 (`IMPLEMENTATION_STATUS.md`, pilot results). We built
two candidates and compared them with a conjugate baseline:

- C0: independent Dirichlet(1/2 + y) per sample, no Stan.
- C1: M1 (`model = "m1"`).
- C2: M1 mean plus a Dirichlet-multinomial sample residual (`model = "m1dm"`).

The run covered 545 replicates and 18 A/B/C cells at 24 donors. Source:
`runs/2026-11-02_redesign-summary/` and `runs/2026-11-02_redesign-decision/`
in the reproducibility repo.

## The pre-registered rule

Committed before any result (outer repo `fed35d9` and the plan):

1. A model passes a cell when RMS, CPS and CCS 90% coverage are each in
   0.85 to 0.97, the false sign rate is at most 0.07, and `ipf_failed` is
   below 1%.
2. A model is calibrated when it passes every A/B/C cell.
3. Among calibrated models, the lowest RMS RMSE on rare_weak plus
   rare_moderate wins.
4. If the best structured model does not beat C0 by 10% on that RMSE, C0 wins.

Amendment 1, before any result (`6e51462`): the band is widened by 2
seed-clustered SEs, and null cells have no upper bound, because every truth
there is exactly 0 and over-coverage is conservative.

Amendment 2, after the results (`de918d2`, post hoc): decide on RMS and CPS
only. CCS becomes experimental. Reason: no candidate covered CCS under
generator B (0.27 to 0.39 across C1 and C2, `calibration.tsv`).

## Result of the rule

With the amended score set the rule still found no calibrated model
(`decision.md`). C0 passed 5 of 18 cells, C1 9 and C2 6. The failures differ in
direction:

| Model | Failing cells | Direction |
|---|---|---|
| C1 (M1) | 9 of 18 | CPS under-covers, 0.24 to 0.79 |
| C2 (m1dm) | 12 of 18 | CPS over-covers, 0.988 to 0.999 (nominal 0.90) |
| C0 | 13 of 18 | mixed; CCS false signs up to 0.99 |

## Decision

`cj_fit()` defaults to `model = "m1dm"`. `"m1"` stays for ablation. The user
chose C2 on 2026-11-02, against the rule's "no default" output.

Why C2: RMS coverage is in the band for C1 and C2 in every cell. C2's
failures are conservative, C1's understate uncertainty. Both cut RMS RMSE on
rare_weak plus rare_moderate by 31% against C0
(`calibration.tsv`).

## Rejected options

- **C1 (m1) as default.** Its CPS intervals are too narrow in 9 cells.
  Understated uncertainty is the worse error.
- **C0 (Dirichlet) as default.** the hierarchies cut RMS RMSE by 31% against it on rare cells, rare-weak
  detection AUC is 0.80 against 0.98 to 0.99 for the hierarchies, and CCS false
  signs reach 0.995 (`auc.tsv`, `false_sign.tsv`).
- **No default until C2's CPS over-coverage is fixed.** Declined: it leaves
  users with no usable model for an unknown wait.

## Consequences

- CPS intervals under `m1dm` are conservative. Report them as such.
- CCS is experimental under every model. Score tables carry an
  `experimental` column and `cj_scores()` says so once per call.
- Fits cost more: about 1.5 times M1's median time (`runtime.tsv`).
- The cause of C2's CPS over-coverage is unknown. The Dirichlet residual may
  inflate ratio-type scores. This is unverified follow-up.
- Held-out generator D was opened only after this decision. Its result is
  recorded in `IMPLEMENTATION_STATUS.md` and can reverse this ADR.
