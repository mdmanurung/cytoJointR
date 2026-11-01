# Visit change: per-draw difference in each net score between two visits of
# the same donor, stimulation and parent population (GLOSSARY: Visit change).

#' Change in net scores between two visits
#'
#' For each donor, stimulation and parent population, subtracts the `from`
#' visit's net score from the `to` visit's, draw by draw. The result is
#' `unavailable` when either visit's contrast has no control sample or the
#' donor lacks a contrast at either visit.
#'
#' The population table averages the per-draw donor changes with equal donor
#' weights over the donors that have both visits available. It describes the
#' observed donors only, not a new donor. A draw where any contributing donor
#' failed is counted as failed.
#'
#' @inheritParams cj_scores
#' @param from,to Visit levels (strings). Change is `to` minus `from`.
#' @return A list with `donor` (one row per donor, stimulation, parent
#'   population, score and marker, with a `reason` for unavailable rows) and
#'   `population` (one row per stimulation, parent population, score and
#'   marker, with `n_donors`, `n_donors_unavailable` and `target` =
#'   `"observed-donor average, equal donor weights"`).
#' @export
#' @examples
#' d <- cj_example_draws(n_draw = 50)
#' eps <- c(RMS = 0.005, poly_absolute = 0.005, poly_marginal = 0.005, CPS = 0.02,
#'          CCS = 0.005, breadth = 0.05, marker_magnitude = 0.005)
#' vc <- cj_visit_change(d, eps, from = "v1", to = "v2")
#' head(vc$population)
cj_visit_change <- function(x, epsilon, from, to, level = 0.9, max_failure = 0.01) {
  d <- as_draws(x)
  check_epsilon(epsilon, c(GLOBAL_SCORES, "marker_magnitude"))
  lv <- levels(d$contrasts$visit)
  if (!all(c(from, to) %in% lv) || from == to)
    stop("from and to must be two different visit levels: ", paste(lv, collapse = ", "),
         call. = FALSE)
  eng <- score_engine(d)
  ct <- d$contrasts
  key <- paste(ct$donor_id, ct$stimulation, ct$parent_population, sep = "\r")
  cols <- eng[[1]]$cols
  eps <- epsilon[split_score_col(cols)$score]
  keys <- unique(data.frame(donor_id = ct$donor_id, stimulation = ct$stimulation,
                            parent_population = ct$parent_population, key = key))
  # Per donor: change matrix and failure flag, or NULL with the reason.
  per_donor <- lapply(keys$key, function(k) {
    i_from <- which(key == k & ct$visit == from)
    i_to <- which(key == k & ct$visit == to)
    if (!length(i_from) || !length(i_to)) return(list(reason = "no contrast at one visit"))
    a <- eng[[i_from]]; b <- eng[[i_to]]
    if (!a$available || !b$available)
      return(list(reason = "no control sample at one visit"))
    list(delta = b$G - a$G, bad = a$status != "ok" | b$status != "ok",
         ids = c(ct$contrast_id[i_from], ct$contrast_id[i_to]))
  })
  donor <- do.call(rbind, lapply(seq_len(nrow(keys)), function(r) {
    pd <- per_donor[[r]]
    ok <- !is.null(pd$delta)
    cbind(keys[r, c("donor_id", "stimulation", "parent_population")],
          from = from, to = to,
          contrast_id_from = if (ok) pd$ids[1] else NA_character_,
          contrast_id_to = if (ok) pd$ids[2] else NA_character_,
          summarise_columns(cols, pd$delta, pd$bad, eps, level, max_failure),
          reason = if (ok) NA_character_ else pd$reason, row.names = NULL)
  }))
  # Population: equal-weight mean over donors that have both visits.
  grp <- paste(keys$stimulation, keys$parent_population, sep = "\r")
  population <- do.call(rbind, lapply(unique(grp), function(g) {
    r <- which(grp == g)
    have <- Filter(function(i) !is.null(per_donor[[i]]$delta), r)
    delta <- if (length(have))
      Reduce(`+`, lapply(have, function(i) per_donor[[i]]$delta)) / length(have)
    bad <- if (length(have)) Reduce(`|`, lapply(have, function(i) per_donor[[i]]$bad))
    cbind(keys[r[1], c("stimulation", "parent_population")], from = from, to = to,
          summarise_columns(cols, delta, bad, eps, level, max_failure),
          n_donors = length(have), n_donors_unavailable = length(r) - length(have),
          target = "observed-donor average, equal donor weights", row.names = NULL)
  }))
  list(donor = donor, population = population)
}
