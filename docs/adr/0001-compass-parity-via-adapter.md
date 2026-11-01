# 0001. COMPASS parity through an adapter plus our own model

Date: 2026-11-01. Status: accepted.

## Context

cytoJointR must offer every feature of COMPASS. COMPASS rests on a responder
mixture model: a latent 0/1 response indicator per subject x subset, with
FS/PFS derived from it. The design bundle lists responder-mixture
identification as a v1 non-goal and plans a donor-paired log-linear model (M1)
with signed contrasts instead.

## Decision

Native COMPASS outputs (FS, PFS, response probabilities, heatmaps) come from
the real `COMPASS::COMPASS()` through an adapter; COMPASS is a `Suggests`
dependency. cytoJointR's own model provides a signed counterpart for each
COMPASS output. We do not reimplement the COMPASS mixture model.

## Consequences

- Native FS/PFS stay exactly COMPASS's numbers and are never relabelled as
  cytoJointR sign probabilities.
- One inference model to validate inside the package, not two.
- The adapter needs COMPASS installed; without it those outputs error clearly.

Rejected: reimplementing COMPASS in Stan (two models, double validation).
