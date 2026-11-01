# Descriptive layer: counts, percentages and signed background subtraction per
# contrast. Written from the spec; clipping exists only in net_percent_display.

#' @importFrom ggplot2 .data
#' @noRd
NULL

#' Describe a `cj_data` object per contrast and pattern
#'
#' @param x A `cj_data` object.
#' @return A `cj_description` list with three data frames:
#' * `patterns`: one row per contrast x pattern with `n_stim`, `pct_stim`,
#'   `n_ctrl`, `pct_ctrl` (percent of the parent count), `net_percent`
#'   (stimulated minus control, signed and never clipped), `net_percent_display`
#'   (`pmax(net_percent, 0)`), `display_only` (`TRUE` where the display value
#'   differs from `net_percent`), `degree` and `functionality` (`none`, `mono`,
#'   `poly`). Contrasts without a control sample have `NA` control and net
#'   columns; samples with `N = 0` have `NA` percentages.
#' * `exactly`: per contrast x degree d, the sums over patterns of degree d.
#' * `at_least`: per contrast x d >= 1, the sums over patterns of degree >= d.
#'
#' The `display_only` value is for plots only and is not used for inference.
#' @examples
#' d <- cj_describe(cj_example())
#' head(d$patterns)
#' d$at_least
#' @export
cj_describe <- function(x) {
  stopifnot(inherits(x, "cj_data"))
  ct <- x$contrasts
  K <- nrow(x$bits)
  si <- match(ct$stim_sample_id, x$samples$sample_id)
  ci <- match(ct$control_sample_id, x$samples$sample_id)  # NA when no control
  rep_k <- function(v) rep(v, each = K)
  by_contrast <- function(m, i) as.vector(t(m[i, , drop = FALSE]))
  pct <- x$counts / x$samples$N * 100  # N = 0 gives NaN
  pct[is.nan(pct)] <- NA
  deg <- rowSums(x$bits)
  p <- data.frame(
    contrast_id = rep_k(ct$contrast_id), donor_id = rep_k(ct$donor_id), visit = rep_k(ct$visit),
    parent_population = rep_k(ct$parent_population), stimulation = rep_k(ct$stimulation),
    pattern = rep(colnames(x$counts), nrow(ct)), degree = rep(as.integer(deg), nrow(ct)),
    n_stim = by_contrast(x$counts, si), pct_stim = by_contrast(pct, si),
    n_ctrl = by_contrast(x$counts, ci), pct_ctrl = by_contrast(pct, ci))
  p$functionality <- c("none", "mono", "poly")[pmin(p$degree, 2L) + 1L]
  p$net_percent <- p$pct_stim - p$pct_ctrl
  p$net_percent_display <- pmax(p$net_percent, 0)
  p$display_only <- p$net_percent_display != p$net_percent & !is.na(p$net_percent)
  p$display_only[is.na(p$net_percent)] <- FALSE

  cols <- c("n_stim", "pct_stim", "n_ctrl", "pct_ctrl", "net_percent")
  agg <- function(keep) {  # keep(degree, d) selects the patterns summed into row d
    out <- do.call(rbind, lapply(sort(unique(p$degree)), function(k) {
      q <- p[keep(p$degree, k), ]
      s <- rowsum(q[cols], q$contrast_id, reorder = FALSE)  # NA stays NA
      data.frame(contrast_id = rownames(s), degree = k, s, row.names = NULL)
    }))
    meta <- ct[match(out$contrast_id, ct$contrast_id), c("donor_id", "visit", "parent_population", "stimulation")]
    out <- cbind(out[1L], meta, out[-1L])
    out[order(match(out$contrast_id, ct$contrast_id), out$degree), ]
  }
  exactly <- agg(`==`)
  at_least <- agg(`>=`)
  at_least <- at_least[at_least$degree >= 1L, ]
  rownames(exactly) <- rownames(at_least) <- NULL
  structure(list(patterns = p, exactly = exactly, at_least = at_least), class = "cj_description")
}

#' Small example `cj_data`
#'
#' Two donors, one visit, markers IFNG and TNF, one stimulation (CMV), each with a control sample.
#' @return A `cj_data` object.
#' @examples
#' cj_example()
#' @export
cj_example <- function() {
  s <- expand.grid(sample_id = c("D1_ctrl", "D1_cmv", "D2_ctrl", "D2_cmv"), IFNG = 0:1, TNF = 0:1,
                   stringsAsFactors = FALSE)
  s$count <- c(900, 800, 950, 700, 50, 90, 30, 120, 30, 60, 15, 100, 20, 50, 5, 80)
  samples <- data.frame(
    sample_id = c("D1_ctrl", "D1_cmv", "D2_ctrl", "D2_cmv"), donor_id = rep(c("D1", "D2"), each = 2),
    visit = factor("v1"), parent_population = "CD4", condition = c("control", "CMV"))
  cj_from_counts(s, samples, c("IFNG", "TNF"))
}

#' Pattern-by-marker heatmap with full state labels
#'
#' One row per exact pattern (labelled with every marker's state, such as
#' `IFNG+ TNF-`), one column per marker, ordered by degree. With `contrast`,
#' a last column prints that contrast's signed `net_percent`.
#'
#' @param x A `cj_data` object.
#' @param contrast Optional `contrast_id` whose signed net percent is printed.
#' @return A ggplot2 object. ComplexHeatmap is not used because the result must
#'   be a ggplot2 object.
#' @examples
#' cj_plot_pattern_heatmap(cj_example(), contrast = "D1_v1_CD4_CMV")
#' @export
cj_plot_pattern_heatmap <- function(x, contrast = NULL) {
  stopifnot(inherits(x, "cj_data"))
  lab <- pattern_labels(x$bits)
  ord <- order(rowSums(x$bits), decreasing = TRUE)
  tiles <- data.frame(pattern = factor(rep(lab, ncol(x$bits)), levels = rev(lab[ord])),
                      marker = factor(rep(x$markers, each = length(lab)), levels = x$markers),
                      state = ifelse(as.vector(x$bits) == 1L, "positive", "negative"))
  g <- ggplot2::ggplot(tiles, ggplot2::aes(.data$marker, .data$pattern)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$state), colour = "white") +
    ggplot2::scale_fill_manual(values = c(positive = "grey20", negative = "grey90")) +
    ggplot2::labs(x = NULL, y = NULL, fill = "State")
  if (is.null(contrast)) return(g)
  d <- cj_describe(x)$patterns
  d <- d[d$contrast_id == contrast, ]
  if (!nrow(d)) stop("unknown contrast: ", contrast, call. = FALSE)
  d$net_label <- ifelse(is.na(d$net_percent), "NA", sprintf("%+.2f", d$net_percent))
  g + ggplot2::geom_text(
    data = data.frame(marker = "net %", pattern = factor(d$pattern, levels = levels(tiles$pattern)),
                      label = d$net_label),
    ggplot2::aes(.data$marker, .data$pattern, label = .data$label), inherit.aes = FALSE) +
    ggplot2::scale_x_discrete(limits = c(x$markers, "net %"))
}

#' Bars of net percent by degree
#'
#' Mean signed net percent per stimulation and degree (bars from 0), with one
#' point per contrast. Values are signed; nothing is clipped.
#'
#' @param x A `cj_data` object.
#' @param type `"exactly"` (patterns of degree d) or `"at_least"` (degree d or more).
#' @return A ggplot2 object. Contrasts without a control sample have no net
#'   percent; their number is reported with a message.
#' @examples
#' cj_plot_degree_bars(cj_example(), "at_least")
#' @export
cj_plot_degree_bars <- function(x, type = c("exactly", "at_least")) {
  type <- match.arg(type)
  d <- cj_describe(x)[[type]]
  if (n_na <- sum(is.na(d$net_percent))) message(n_na, " contrast x degree rows have no control sample and are not plotted")
  d <- d[!is.na(d$net_percent), ]
  m <- stats::aggregate(net_percent ~ stimulation + degree, d, mean)
  ggplot2::ggplot(m, ggplot2::aes(factor(.data$degree), .data$net_percent)) +
    ggplot2::geom_col(fill = "grey70") +
    ggplot2::geom_point(data = d, position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1)) +
    ggplot2::facet_wrap(~stimulation) +
    ggplot2::labs(x = if (type == "exactly") "Degree (exactly)" else "Degree (at least)",
                  y = "Net percent of parent (stimulated - control)")
}

#' Two-marker scatters coloured by pattern
#'
#' For every marker pair, plots the cells' fluorescence coloured by their exact
#' pattern. Patterns come from `thresholds` (above is positive).
#'
#' @param cells Per-cell data frame with `sample_id` and one fluorescence column per marker.
#' @param markers Marker names (at least two).
#' @param thresholds Named numeric vector, names equal to `markers`.
#' @param sample_id Optional single sample to plot; `NULL` pools all cells.
#' @param max_cells Cells drawn at most; a larger table is randomly thinned with a message.
#' @param seed Seed for the thinning.
#' @return A ggplot2 object.
#' @examples
#' cells <- data.frame(sample_id = "s", A = c(0, 2, 2, 0, 3), B = c(0, 0, 2, 3, 3))
#' cj_plot_pattern_scatter(cells, c("A", "B"), c(A = 1, B = 1))
#' @export
cj_plot_pattern_scatter <- function(cells, markers, thresholds, sample_id = NULL, max_cells = 20000L,
                                    seed = 1L) {
  stopifnot(is.data.frame(cells), length(markers) >= 2L, identical(names(thresholds), markers),
            all(markers %in% names(cells)))
  if (!is.null(sample_id)) {
    cells <- cells[cells$sample_id == sample_id, ]
    if (!nrow(cells)) stop("no cells for sample_id ", sample_id, call. = FALSE)
  }
  if (nrow(cells) > max_cells) {
    message("thinning ", nrow(cells), " cells to ", max_cells, " for plotting (seed ", seed, ")")
    keep <- withr_seed(seed, sample.int(nrow(cells), max_cells))
    cells <- cells[keep, ]
  }
  z <- sweep(as.matrix(cells[markers]), 2L, thresholds, `>`)
  pattern <- apply(z, 1L, function(b) paste0(markers, ifelse(b, "+", "-"), collapse = " "))
  pairs <- utils::combn(markers, 2L, simplify = FALSE)
  d <- do.call(rbind, lapply(pairs, function(p)
    data.frame(x = cells[[p[1]]], y = cells[[p[2]]], pair = paste(p, collapse = " vs "), pattern = pattern)))
  ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y, colour = .data$pattern)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.5) +
    ggplot2::facet_wrap(~pair, scales = "free") +
    ggplot2::labs(x = NULL, y = NULL, colour = "Pattern")
}

# Run expr with a seed without changing the caller's RNG stream.
withr_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, globalenv()))
  set.seed(seed)
  expr
}

#' Print and summarise a `cj_data` object
#'
#' @param x,object A `cj_data` object.
#' @param ... Unused.
#' @return `print` returns `x` invisibly; `summary` returns a `summary.cj_data`
#'   list (donors, visits, stimulations, contrasts without a control sample, N
#'   range) that prints.
#' @examples
#' x <- cj_example()
#' x
#' summary(x)
#' @export
summary.cj_data <- function(object, ...) {
  ct <- object$contrasts
  structure(list(
    donors = unique(object$samples$donor_id), visits = levels(object$samples$visit),
    stimulations = unique(ct$stimulation), n_samples = nrow(object$samples),
    n_contrasts = nrow(ct), no_control = ct$contrast_id[is.na(ct$control_sample_id)],
    N_range = range(object$samples$N), markers = object$markers), class = "summary.cj_data")
}

#' @rdname summary.cj_data
#' @export
print.summary.cj_data <- function(x, ...) {
  cat("cj_data: ", x$n_samples, " samples, ", length(x$markers), " markers (", 2^length(x$markers),
      " patterns)\n", sep = "")
  cat("Donors (", length(x$donors), "): ", paste(x$donors, collapse = ", "), "\n", sep = "")
  cat("Visits: ", paste(x$visits, collapse = ", "), "\n", sep = "")
  cat("Stimulations: ", paste(x$stimulations, collapse = ", "), "\n", sep = "")
  cat("Contrasts: ", x$n_contrasts, "; without a control sample: ", length(x$no_control),
      if (length(x$no_control)) paste0(" (", paste(x$no_control, collapse = ", "), ")"), "\n", sep = "")
  cat("N range: ", x$N_range[1], " to ", x$N_range[2], "\n", sep = "")
  invisible(x)
}

#' @rdname summary.cj_data
#' @export
print.cj_data <- function(x, ...) {
  print(summary(x))
  invisible(x)
}
