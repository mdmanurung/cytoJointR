# cytoJointR

Donor-paired inference on polyfunctional intracellular cytokine/function assays:
how stimulation changes the distribution of exact marker patterns within donors,
relative to each donor's matched control sample.

## Samples and design

**Donor**:
The biological replicate; one person or animal contributing samples.
_Avoid_: subject, individual, patient (COMPASS's "subject" maps here)

**Visit**:
One sampling occasion of a donor; all samples of a donor at one visit share
that visit's control sample.
_Avoid_: timepoint, day

**Parent population**:
The gated cell population whose cells are classified into patterns, such as CD4+ T cells.
_Avoid_: subset, cell subset, cell type

**Stimulation**:
A named treatment applied to an aliquot, such as CMV peptide pool or SEB.
_Avoid_: treatment, antigen, condition

**Control sample**:
The unstimulated aliquot of a donor at a visit; it may serve several stimulations.
_Avoid_: unstim, background sample, negative control

**Contrast**:
One stimulated sample paired with its control sample; the unit a net score describes.
_Avoid_: pair, comparison

## Patterns and markers

**Functional marker**:
A measured cell function such as IFNG, TNF, IL2 or CD107a.
_Avoid_: cytokine (CD107a is not one)

**Pattern**:
One exact Boolean combination over all functional markers, including negative states;
every cell falls in exactly one.
_Avoid_: subset, combination, category, Boolean gate

**Degree**:
The number of positive markers in a pattern.
_Avoid_: functionality, order

**Polyfunctional**:
A pattern with degree 2 or more.

## Probabilities and scores

**Background**:
The fitted pattern probabilities of a control sample.
_Avoid_: baseline, reference, unstimulated frequency

**Net score**:
A signed stimulated-minus-background quantity for one contrast, never clipped.
_Avoid_: background-subtracted frequency, response

**Visit change**:
The difference in one net score between two visits of the same donor and stimulation.
_Avoid_: longitudinal change, delta-delta

**Reference**:
A fixed comparison distribution used only to define a score: the counterfactual
for CCS or the fixed pattern weights for FIS.
_Avoid_: baseline, control

**Display-only value**:
A derived value, such as a clipped net score, shown in plots but never used for inference.
