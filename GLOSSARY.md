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

**Missing control**:
A stimulated sample with no control sample at its visit. Net scores stay
unavailable. COMPASS can still run it (`cj_compass(missing_control = )`): with
the same donor's control from the nearest other visit (`"donor"`), or against
the population control of the other donors (`"population"`); `control_source`
says which.
_Avoid_: unpaired sample, orphan

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

## Fluorescence and gating

**Pooled gate**:
One threshold per marker, found on all control events pooled across donors,
visits and samples, on a stated transform (`cj_gate()`).
_Avoid_: per-sample gate, cutoff

**Threshold** (`t_j`):
The transformed-fluorescence value above which marker `j` is positive.

**Excess** (`X_j - t_j`):
Transformed fluorescence above the threshold. Positive for a positive cell.
Units are those of the transform (for example asinh units), not protein counts.

**Conditional intensity** (`I_ckj`):
The mean excess of marker `j` among the cells of pattern `k` in condition `c`.
Plug-in estimate from events; no truncated-Gaussian model in the descriptive scores.

**Burden** (`B_cj`):
`sum over patterns positive for j of p_ck I_ckj`: the mean of `(X_j - t_j) * 1[X_j > t_j]`
over all parent events of condition `c`.

**NICR** (net integrated cytokine response, `NICR_j = B_Sj - B_Uj`):
Signed stimulated-minus-control burden of one marker. Never summed across
markers, because channels are on different scales.
_Avoid_: iMFI (a different, mean- or median-based quantity)

**Poly-NICR**:
The NICR sum restricted to patterns of degree 2 or more.

**Frequency and intensity parts**:
Exact split `dB_k = dp_k (I_S + I_U)/2 + dI_k (p_S + p_U)/2` per pattern. A
pattern seen in one condition only goes to `nicr_unmatched`. An algebraic
allocation, not a causal one.

**FIS** (fixed-reference intensity shift):
`sum_k w_kj (I_Skj - I_Ukj)` with weights `w_kj` fixed and common to both
conditions. Unavailable (`NA` with a reason) where a reference pattern lacks events.

**Intensity-weighted FS and PFS** (`iw_fs`, `iw_pfs`; experimental):
COMPASS's FS and PFS weights applied to the net burden of marker `j` in each
subset, `dB_kj = p_Sk I_Skj - p_Uk I_Ukj`, in place of COMPASS's posterior
response probability. `iw_fs_j = sum dB_kj / (2^J - 1)` (equal to `NICR_j` over
`2^J - 1`); `iw_pfs_j = sum dB_kj * d_k / choose(J, d_k) / J * 2 / (J + 1)`.
Signed, per marker, in fluorescence units per parent cell; not on the scale of
COMPASS FS/PFS.
