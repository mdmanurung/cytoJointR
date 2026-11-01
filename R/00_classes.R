# Core S3 classes shared by the data, model and score layers.
# Contract (docs/SCOPE.md): every layer passes these objects; nobody invents a
# second class for the same thing. Vocabulary follows GLOSSARY.md.
#
# cj_data   counts [sample x pattern] integer matrix, all 2^J patterns
#           samples   one row per sample, row order == rownames(counts)
#           contrasts one row per stimulated sample; control_sample_id NA when
#                     the visit has no control sample (net scores unavailable)
#           bits      K x J integer matrix from cj_bits(), first marker low bit
# cj_draws  p [draw x sample x pattern] fitted pattern probabilities, plus the
#           same samples / contrasts / bits, plus engine metadata
# cj_fit    data + draws + engine object + design + priors + diagnostics

CONTROL_CONDITION <- "control"

SAMPLE_COLUMNS <- c("sample_id", "donor_id", "visit", "parent_population",
                    "condition", "N")
CONTRAST_COLUMNS <- c("contrast_id", "donor_id", "visit", "parent_population",
                      "stimulation", "stim_sample_id", "control_sample_id")
DRAW_KINDS <- c("in_sample_donor_conditional", "new_donor", "reference")

# Full state label per pattern, e.g. "IFNG+ TNF- IL2-"; shared by data and score layers.
pattern_labels <- function(bits) {
  apply(bits, 1L, function(b)
    paste0(colnames(bits), ifelse(b == 1L, "+", "-"), collapse = " "))
}

new_cj_data <- function(counts, samples, contrasts, bits, provenance = list()) {
  validate_cj_data(structure(
    list(counts = counts, samples = samples, contrasts = contrasts,
         bits = bits, markers = colnames(bits), provenance = provenance),
    class = "cj_data"))
}

validate_cj_data <- function(x) {
  validate_design_parts(x$samples, x$contrasts, x$bits)
  counts <- x$counts
  if (!is.matrix(counts) || !is.integer(counts) || anyNA(counts) || any(counts < 0L))
    stop("counts must be a nonnegative integer matrix without NA", call. = FALSE)
  if (!identical(rownames(counts), x$samples$sample_id))
    stop("rownames(counts) must equal samples$sample_id in order", call. = FALSE)
  if (!identical(colnames(counts), pattern_labels(x$bits)))
    stop("colnames(counts) must be pattern_labels(bits) in pattern order", call. = FALSE)
  if (!identical(as.numeric(rowSums(counts)), as.numeric(x$samples$N)))
    stop("samples$N must equal the row sums of counts", call. = FALSE)
  x
}

new_cj_draws <- function(p, samples, contrasts, bits, meta) {
  validate_cj_draws(structure(
    list(p = p, samples = samples, contrasts = contrasts, bits = bits,
         markers = colnames(bits), meta = meta),
    class = "cj_draws"))
}

validate_cj_draws <- function(x, tol = 1e-8) {
  validate_design_parts(x$samples, x$contrasts, x$bits)
  p <- x$p
  if (!is.array(p) || length(dim(p)) != 3L)
    stop("p must be a [draw, sample, pattern] array", call. = FALSE)
  if (!identical(dimnames(p)[[2]], x$samples$sample_id))
    stop("dimnames(p)[[2]] must equal samples$sample_id in order", call. = FALSE)
  if (!identical(dimnames(p)[[3]], pattern_labels(x$bits)))
    stop("dimnames(p)[[3]] must be pattern_labels(bits)", call. = FALSE)
  if (any(!is.finite(p)) || any(p <= 0))
    stop("p must be finite and strictly positive (model probabilities have full support)",
         call. = FALSE)
  if (max(abs(apply(p, c(1L, 2L), sum) - 1)) > tol)
    stop("p must sum to one over patterns in every draw and sample", call. = FALSE)
  meta_fields <- c("engine_id", "model_id", "seed", "n_draws", "draw_kind")
  if (!all(meta_fields %in% names(x$meta)))
    stop("meta needs: ", paste(meta_fields, collapse = ", "), call. = FALSE)
  if (!x$meta$draw_kind %in% DRAW_KINDS)
    stop("meta$draw_kind must be one of: ", paste(DRAW_KINDS, collapse = ", "),
         call. = FALSE)
  if (x$meta$n_draws != dim(p)[1])
    stop("meta$n_draws must equal dim(p)[1]", call. = FALSE)
  x
}

new_cj_fit <- function(data, draws, engine_fit, design, priors, diagnostics,
                       spec_version = "M1-0.1") {
  stopifnot(inherits(data, "cj_data"), inherits(draws, "cj_draws"))
  if (!identical(data$samples, draws$samples) ||
      !identical(data$contrasts, draws$contrasts) ||
      !identical(data$bits, draws$bits))
    stop("draws must describe the same samples, contrasts and patterns as data",
         call. = FALSE)
  structure(list(data = data, draws = draws, engine_fit = engine_fit,
                 design = design, priors = priors, diagnostics = diagnostics,
                 spec_version = spec_version),
            class = "cj_fit")
}

# Shared checks on the design tables; both cj_data and cj_draws carry them.
validate_design_parts <- function(samples, contrasts, bits) {
  if (!is.matrix(bits) || !is.integer(bits) || is.null(colnames(bits)) ||
      nrow(bits) != 2L^ncol(bits) ||
      !identical(unname(bits), unname(cj_bits(colnames(bits)))))
    stop("bits must be cj_bits(markers)", call. = FALSE)
  missing_s <- setdiff(SAMPLE_COLUMNS, names(samples))
  missing_c <- setdiff(CONTRAST_COLUMNS, names(contrasts))
  if (length(missing_s)) stop("samples lacks: ", paste(missing_s, collapse = ", "), call. = FALSE)
  if (length(missing_c)) stop("contrasts lacks: ", paste(missing_c, collapse = ", "), call. = FALSE)
  if (anyDuplicated(samples$sample_id)) stop("duplicate sample_id", call. = FALSE)
  if (anyDuplicated(contrasts$contrast_id)) stop("duplicate contrast_id", call. = FALSE)
  if (!is.factor(samples$visit) || !is.factor(contrasts$visit))
    stop("visit must be a factor with explicit levels", call. = FALSE)
  if (any(contrasts$stimulation == CONTROL_CONDITION))
    stop("a stimulation may not be named '", CONTROL_CONDITION, "'", call. = FALSE)

  row_of <- match(contrasts$stim_sample_id, samples$sample_id)
  if (anyNA(row_of)) stop("contrasts reference unknown stim_sample_id", call. = FALSE)
  if (!identical(samples$condition[row_of], contrasts$stimulation))
    stop("stim_sample_id condition must equal the contrast's stimulation", call. = FALSE)
  key <- c("donor_id", "parent_population")
  for (k in key)
    if (!identical(as.character(samples[[k]][row_of]), as.character(contrasts[[k]])))
      stop("contrast ", k, " disagrees with its stimulated sample", call. = FALSE)
  if (!identical(as.character(samples$visit[row_of]), as.character(contrasts$visit)))
    stop("contrast visit disagrees with its stimulated sample", call. = FALSE)

  has_ctrl <- !is.na(contrasts$control_sample_id)
  ctrl_row <- match(contrasts$control_sample_id[has_ctrl], samples$sample_id)
  if (anyNA(ctrl_row)) stop("contrasts reference unknown control_sample_id", call. = FALSE)
  if (any(samples$condition[ctrl_row] != CONTROL_CONDITION))
    stop("control_sample_id must point at a control sample", call. = FALSE)
  for (k in key)
    if (!identical(as.character(samples[[k]][ctrl_row]),
                   as.character(contrasts[[k]][has_ctrl])))
      stop("control sample ", k, " disagrees with its contrast", call. = FALSE)
  if (!identical(as.character(samples$visit[ctrl_row]),
                 as.character(contrasts$visit[has_ctrl])))
    stop("control sample must come from the contrast's visit", call. = FALSE)
  invisible(TRUE)
}
