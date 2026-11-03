# One-call wrappers for the adoption API: from a cj_data (or straight from cells)
# to tidy per-donor tables. They reuse cj_compass() and the polyICSFlow
# functions; optional packages are checked with has_pkg().

# A cj_data from either a cj_data or a cell data frame.
as_cj_data <- function(x, markers, samples, gate, transform, control, ...) {
  if (inherits(x, "cj_fit")) return(x$data)
  if (inherits(x, "cj_data")) return(x)
  if (!is.data.frame(x)) stop("x must be a cj_data or a cell data frame", call. = FALSE)
  if (is.null(samples)) stop("a cell data frame needs samples", call. = FALSE)
  cj_from_cells(x, markers, samples, gate = gate, transform = transform, control = control, ...)
}

#' Run COMPASS in one call and return tidy FS and PFS per donor
#'
#' A thin layer over [cj_compass()]. Starting from cells, it first runs
#' [cj_from_cells()] (and [cj_gate()] if no gate is given). The table has the
#' COMPASS quantities `compass_FS`, `compass_PFS` and `compass_response_degree1`,
#' `compass_response_degree2`, defined in [cj_compass()]. They are COMPASS
#' numbers and not interchangeable with cytoJointR sign probabilities.
#'
#' @param x A `cj_data`, a `cj_fit`, or a cell data frame.
#' @param markers,samples,gate,transform,control Used only when `x` is a cell
#'   data frame; see [cj_from_cells()]. `samples` is required then.
#' @param missing_control What to do with stimulated samples that have no
#'   control: `"drop"` (default), `"donor"` or `"population"`; see [cj_compass()].
#'   The table gets `control_source` and `control_fallback`.
#' @param ... Passed to [cj_compass()] (for example `iterations`, `replications`,
#'   `seed`).
#' @param gate_args List of extra arguments for [cj_gate()] (for example
#'   `list(method = "quantile")`) when `x` is a cell data frame without a `gate`.
#' @return A data frame, one row per donor, stimulation, visit and parent
#'   population, with the attributes `compass` (the full `cj_compass` object)
#'   and `excluded` (contrasts left out, with reason).
#' @export
#' @examplesIf requireNamespace("COMPASS", quietly = TRUE)
#' \donttest{
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 400, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' samples$visit <- factor("v1")
#' samples$parent_population <- "CD4"
#' cj_run_compass(sim$cells, c("IFNG", "TNF", "IL2"), samples, transform = "none",
#'                control = "unstim", gate_args = list(method = "quantile"),
#'                iterations = 200, replications = 2)
#' }
cj_run_compass <- function(x, markers = NULL, samples = NULL, gate = NULL, transform = NULL,
                           control = "control", missing_control = c("drop", "donor", "population"),
                           ..., gate_args = list()) {
  if (!has_pkg("COMPASS")) stop("cj_run_compass() needs the COMPASS package", call. = FALSE)
  d <- do.call(as_cj_data, c(list(x, markers, samples, gate, transform, control), gate_args))
  cc <- cj_compass(d, missing_control = match.arg(missing_control), ...)
  out <- cc$scores
  attr(out, "compass") <- cc
  attr(out, "excluded") <- cc$excluded
  if (nrow(cc$excluded))
    message(nrow(cc$excluded), " contrast(s) were not run through COMPASS (see attr(, \"excluded\") for the reason)")
  out
}

# polyICSFlow's own marker-combination labels, e.g. "IFNG+TNF-".
polyics_labels <- function(bits) {
  apply(bits, 1L, function(b) paste0(colnames(bits), ifelse(b == 1L, "+", "-"), collapse = ""))
}

# One row per event, as polyICSFlow::assignMarkerCombinations() wants it,
# rebuilt from the pattern counts so it agrees with the cj_data exactly.
polyics_events <- function(x) {
  cnt <- as.vector(x$counts)
  s <- rep(row(x$counts)[seq_along(cnt)], cnt)
  k <- rep(col(x$counts)[seq_along(cnt)], cnt)
  pos <- as.data.frame(x$bits[k, , drop = FALSE] == 1L)
  names(pos) <- x$markers
  list(positivity = pos, sample = s)
}

#' Run polyICSFlow in one call and return its net tables
#'
#' Rebuilds one row per event from the pattern counts of a `cj_data` (so the
#' counts are exactly those the other scores use), runs
#' `polyICSFlow::assignMarkerCombinations()` and
#' `polyICSFlow::calcPolyfunctionality()`, and sums the background-subtracted
#' percentages. Needs the polyICSFlow package.
#'
#' polyICSFlow clips each pattern's background-subtracted percentage at 0
#' (`pmax`). `polyics_net_any` and `polyics_net_poly` sum those clipped values
#' over patterns with at least 1 or 2 markers. They are polyICSFlow's quantities.
#' For signed net percentages use [cj_describe()].
#'
#' @inheritParams cj_run_compass
#' @param ... Passed to `polyICSFlow::calcPolyfunctionality()`.
#' @return A list of class `cj_polyicsflow`: `patterns` (polyICSFlow's table per
#'   donor, visit, population, stimulation and marker combination, with
#'   `degree`) and `net` (one row per donor, visit, population and stimulation:
#'   `polyics_net_any`, `polyics_net_poly`, in percent of the parent events).
#'   Contrasts without a control sample have `NA` in `net`.
#' @export
#' @examplesIf requireNamespace("polyICSFlow", quietly = TRUE)
#' \donttest{
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 400, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' samples$visit <- factor("v1")
#' samples$parent_population <- "CD4"
#' cj_run_polyicsflow(sim$cells, c("IFNG", "TNF", "IL2"), samples, transform = "none",
#'                    control = "unstim", gate_args = list(method = "quantile"))$net
#' }
cj_run_polyicsflow <- function(x, markers = NULL, samples = NULL, gate = NULL, transform = NULL,
                               control = "control", ..., gate_args = list()) {
  if (!has_pkg("polyICSFlow")) stop("cj_run_polyicsflow() needs the polyICSFlow package", call. = FALSE)
  d <- do.call(as_cj_data, c(list(x, markers, samples, gate, transform, control), gate_args))
  ev <- polyics_events(d)
  assigned <- polyICSFlow::assignMarkerCombinations(ev$positivity, mode = "simple")
  sm <- d$samples[ev$sample, ]
  cells <- data.frame(donor_id = sm$donor_id, visit = sm$visit, parent_population = sm$parent_population,
                      Stimulation = sm$condition, assigned[c("MarkerComb", "MarkerCount", "Functionality")])
  pat <- polyICSFlow::calcPolyfunctionality(cells, md_cols = c("donor_id", "visit"),
                                            pop_col = "parent_population", resolution = "MarkerComb",
                                            condition_col = "Stimulation",
                                            background_val = CONTROL_CONDITION, ...)
  pat <- as.data.frame(pat)
  pat$degree <- rowSums(d$bits)[match(as.character(pat$MarkerComb), polyics_labels(d$bits))]
  st <- pat[pat$Stimulation != CONTROL_CONDITION, , drop = FALSE]
  net <- function(min_degree) {
    a <- stats::aggregate(perc_pos_backgroundSubtracted ~ donor_id + visit + parent_population + Stimulation,
                          data = st[st$degree >= min_degree, ], FUN = sum, na.action = stats::na.pass)
    a
  }
  any1 <- net(1L); any2 <- net(2L)
  out <- merge(any1, any2, by = c("donor_id", "visit", "parent_population", "Stimulation"),
               suffixes = c("_any", "_poly"), all.x = TRUE)
  names(out) <- c("donor_id", "visit", "parent_population", "stimulation", "polyics_net_any", "polyics_net_poly")
  structure(list(patterns = pat, net = out), class = "cj_polyicsflow")
}

#' Cells to every score in one call
#'
#' A thin pipeline: [cj_gate()] (or the gate you give), [cj_from_cells()],
#' [cj_describe()], [cj_fluorescence_scores()], then COMPASS and polyICSFlow
#' when installed. The Bayesian cytoJointR scores need a Stan fit and are run
#' only with `fit = TRUE`.
#'
#' @inheritParams cj_from_cells
#' @param markers Marker names mapped to columns, as in [cj_gate()].
#' @param n_boot Bootstrap replicates for [cj_fluorescence_scores()].
#' @param run Which optional runs to attempt: `"compass"`, `"polyicsflow"`. A
#'   missing package gives a message, not an error.
#' @param fit Also fit the count model with [cj_fit()] and summarise it with
#'   [cj_scores()]; needs CmdStan and takes minutes.
#' @param epsilon Score thresholds for [cj_scores()], required when `fit = TRUE`.
#' @param fit_args List of further arguments to [cj_fit()].
#' @param compass_args List of arguments to [cj_compass()].
#' @return A list: `gate`, `data`, `description`, `fluorescence`, and, when run,
#'   `compass`, `polyicsflow`, `fit`, `scores`.
#' @export
#' @examples
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 300, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' samples$visit <- factor("v1")
#' samples$parent_population <- "CD4"
#' qs <- cj_quickstart(sim$cells, c("IFNG", "TNF", "IL2"), samples, transform = "none",
#'                     control = "unstim", run = character(), method = "quantile")
#' head(qs$fluorescence[c("donor_id", "marker", "nicr")])
cj_quickstart <- function(cells, markers, samples, gate = NULL, transform = NULL, control = "control",
                          n_boot = 0L, run = c("compass", "polyicsflow"), fit = FALSE, epsilon = NULL,
                          fit_args = list(), compass_args = list(), ...) {
  x <- cj_from_cells(cells, markers, samples, gate = gate, transform = transform, control = control, ...)
  out <- list(gate = x$gate, data = x, description = cj_describe(x),
              fluorescence = cj_fluorescence_scores(x, n_boot = n_boot))
  for (r in intersect(run, c("compass", "polyicsflow"))) {
    pkg <- c(compass = "COMPASS", polyicsflow = "polyICSFlow")[[r]]
    if (!has_pkg(pkg)) { message(pkg, " is not installed; skipping ", r); next }
    out[[r]] <- if (r == "compass") do.call(cj_run_compass, c(list(x), compass_args)) else cj_run_polyicsflow(x)
  }
  if (fit) {
    if (is.null(epsilon)) stop("epsilon is required when fit = TRUE", call. = FALSE)
    out$fit <- do.call(cj_fit, c(list(x), fit_args))
    out$scores <- cj_scores(out$fit, epsilon)
  }
  out
}
