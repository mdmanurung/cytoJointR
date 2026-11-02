#' @importFrom ggplot2 .data
NULL

# Plots. Signed values are never clipped: axes and colour scales are symmetric
# around zero and cover the full range of the data.

# Order of patterns for display: by degree, then by label.
pattern_order <- function(p) {
  u <- unique(p[, c("pattern", "degree")])
  u$pattern[order(u$degree, u$pattern)]
}

#' Signed forest plot of scores per contrast
#'
#' One panel per score (and marker). Points are posterior means, bars are
#' credible intervals, dashed lines are `+-epsilon`. Unavailable contrasts
#' are not drawn and are counted in the caption.
#'
#' @param scores Output of [cj_scores()] or [cj_markers()].
#' @param score Scores to show.
#' @return A ggplot object.
#' @export
#' @examples
#' d <- cj_example_draws(n_draw = 50)
#' eps <- c(RMS = 0.005, poly_absolute = 0.005, poly_marginal = 0.005, CPS = 0.02,
#'          CCS = 0.005, breadth = 0.05, marker_magnitude = 0.005)
#' cj_plot_forest(cj_scores(d, eps), score = c("RMS", "CCS"))
cj_plot_forest <- function(scores, score = c("RMS", "poly_absolute", "CPS", "CCS")) {
  s <- scores[scores$score %in% score, ]
  s$panel <- ifelse(is.na(s$marker), s$score, paste0(s$score, ": ", s$marker))
  n_na <- length(unique(s$contrast_id[s$status == "unavailable"]))
  s <- s[s$status != "unavailable" & !is.na(s$estimate), ]
  ggplot2::ggplot(s, ggplot2::aes(x = .data$estimate, y = .data$contrast_id)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = .data$epsilon), linetype = "dashed",
                        colour = "grey60") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = -.data$epsilon), linetype = "dashed",
                        colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                           orientation = "y", width = 0.25) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ panel, scales = "free_x") +
    ggplot2::labs(x = "net score (posterior mean and credible interval)", y = NULL,
                  caption = paste0("Unavailable contrasts (no control sample), not shown: ", n_na))
}

# Shared prep for the two pattern plots: one contrast per column label.
pattern_plot_data <- function(patterns, score) {
  if (length(unique(patterns$score)) > 1L && is.null(score))
    stop("patterns holds several scores; set score", call. = FALSE)
  if (!is.null(score)) patterns <- patterns[patterns$score == score, ]
  patterns <- patterns[patterns$status != "unavailable", ]
  patterns$pattern <- factor(patterns$pattern, levels = rev(pattern_order(patterns)))
  patterns
}

#' Donor by pattern heatmap
#'
#' Rows are exact patterns with full state labels, ordered by degree (highest
#' degree at the top). Columns are donors, split by stimulation and visit.
#' This is the counterpart of the COMPASS mean-gamma heatmap, but shows
#' signed posterior quantities: the probability that the pattern's
#' contribution exceeds epsilon, or its posterior mean.
#'
#' @param patterns Output of [cj_patterns()].
#' @param value `"p_positive"` for P(contribution > epsilon) or `"estimate"`
#'   for the posterior mean contribution.
#' @param score Which score of `patterns` to draw.
#' @return A ggplot object. The mean scale is symmetric around zero.
#' @export
#' @examples
#' pt <- cj_patterns(cj_example_draws(n_draw = 50), "RMS", c(RMS = 0.005))
#' cj_plot_response_heatmap(pt, "estimate")
cj_plot_response_heatmap <- function(patterns, value = c("p_positive", "estimate"),
                                    score = NULL) {
  value <- match.arg(value)
  pd <- pattern_plot_data(patterns, score)
  g <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$donor_id, y = .data$pattern,
                                        fill = .data[[value]])) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(~ stimulation + visit, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = "donor", y = NULL)
  if (value == "p_positive")
    return(g + ggplot2::scale_fill_viridis_c(limits = c(0, 1), name = "P(contribution > epsilon)"))
  m <- max(abs(pd$estimate), na.rm = TRUE)
  g + ggplot2::scale_fill_gradient2(limits = c(-m, m), name = "posterior mean")
}

#' Diverging pattern-contribution bars
#'
#' Posterior mean contribution of each exact pattern with its credible
#' interval, one panel per contrast. Bars go left for negative and right for
#' positive contributions. Intervals do not add up to the global interval.
#'
#' @inheritParams cj_plot_response_heatmap
#' @param contrast Contrast ids to show; default all.
#' @return A ggplot object.
#' @export
#' @examples
#' pt <- cj_patterns(cj_example_draws(n_draw = 50), "RMS", c(RMS = 0.005))
#' cj_plot_contributions(pt, contrast = "D1_v1_CD4_CMV")
cj_plot_contributions <- function(patterns, score = NULL, contrast = NULL) {
  pd <- pattern_plot_data(patterns, score)
  if (!is.null(contrast)) pd <- pd[pd$contrast_id %in% contrast, ]
  pd$direction <- ifelse(pd$estimate >= 0, "positive", "negative")
  ggplot2::ggplot(pd, ggplot2::aes(x = .data$estimate, y = .data$pattern,
                                   fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                           orientation = "y", width = 0.3) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::scale_fill_manual(values = c(negative = "#b2182b", positive = "#2166ac")) +
    ggplot2::facet_wrap(~ contrast_id) +
    ggplot2::labs(x = "signed contribution to the global score", y = NULL, fill = NULL)
}
