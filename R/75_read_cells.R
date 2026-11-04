# Readers that turn FCS files or a GatingSet node into the per-event `cells` table
# that cj_gate(), cj_from_cells() and cj_quickstart() take.

# Channel of each marker in one flowFrame or cytoframe: `markers` values match a
# channel name or a channel description ($PnS) exactly. Names are the output columns.
fcs_channels <- function(fr, markers) {
  par <- flowCore::parameters(fr)@data
  out <- vapply(markers, function(m) {
    hit <- unique(c(par$name[par$name %in% m], par$name[par$desc %in% m]))
    if (length(hit) != 1L)
      stop("marker '", m, "' matches ", length(hit), " channels; give a channel name or an exact description. ",
           "Channels: ", paste0(par$name, ifelse(is.na(par$desc), "", paste0(" (", par$desc, ")")), collapse = ", "),
           call. = FALSE)
    hit
  }, character(1))
  names(out) <- if (is.null(names(markers))) markers else names(markers)
  out
}

marker_map <- function(markers) {
  if (!is.character(markers) || !length(markers) || anyNA(markers))
    stop("markers must be a character vector of channel names or descriptions", call. = FALSE)
  if (is.null(names(markers))) names(markers) <- markers
  if (any(names(markers) == "") || anyDuplicated(names(markers)))
    stop("marker names must be unique and non-empty", call. = FALSE)
  markers
}

#' Read pre-gated FCS files into a cells table
#'
#' Reads one FCS file per sample with `flowCore::read.FCS()`, applies the
#' file's own spillover matrix and returns the events of the chosen markers in
#' the long `cells` format that [cj_gate()], [cj_from_cells()] and
#' [cj_quickstart()] take.
#'
#' **The files must hold the parent population only** (for example CD4 T cells
#' exported from FlowJo), because `cj_from_cells()` treats every event as a
#' parent event. For ungated files, gate the parent first (for example with an
#' openCyto template on a GatingSet) and use [cj_cells_from_gatingset()].
#'
#' Values are compensated but not transformed; choose the transform in
#' [cj_gate()] (for example `transform = "asinh", cofactor = "auto"`).
#'
#' @param files Paths of the FCS files, one per sample.
#' @param markers Character vector of channel names (`$PnN`, e.g. `"B710-A"`) or
#'   exact channel descriptions (`$PnS`, e.g. `"IFNG BV421"`). Names, when given,
#'   become the output columns (e.g. `c(IFNG = "IFNG BV421")`).
#' @param sample_id Sample identifiers, one per file. Defaults to the file name
#'   without `.fcs`.
#' @param compensate Apply the spillover matrix stored in each file (keyword
#'   `SPILL`, `$SPILLOVER` or `SPILLOVER`, the first one present). A file
#'   without one is an error; pass `FALSE` for files that are already
#'   compensated.
#' @return A data frame with `sample_id` and one numeric column per marker,
#'   with attribute `compensation`: one row per file with the keyword used
#'   (`NA` when `compensate = FALSE`).
#' @export
#' @examplesIf requireNamespace("flowCore", quietly = TRUE) && requireNamespace("flowWorkspaceData", quietly = TRUE)
#' f <- system.file("extdata", c("CytoTrol_CytoTrol_1.fcs", "CytoTrol_CytoTrol_2.fcs"),
#'                  package = "flowWorkspaceData")
#' cells <- cj_read_fcs(f, c(CD38 = "CD38 APC", HLADR = "HLA-DR V500"))
#' head(cells)
#' attr(cells, "compensation")
cj_read_fcs <- function(files, markers, sample_id = NULL, compensate = TRUE) {
  if (!requireNamespace("flowCore", quietly = TRUE))
    stop("cj_read_fcs() needs the flowCore package", call. = FALSE)
  markers <- marker_map(markers)
  if (is.null(sample_id)) sample_id <- sub("\\.fcs$", "", basename(files), ignore.case = TRUE)
  if (length(sample_id) != length(files) || anyDuplicated(sample_id))
    stop("sample_id must give one unique id per file", call. = FALSE)
  if (length(miss <- files[!file.exists(files)])) stop("files not found: ", paste(miss, collapse = ", "), call. = FALSE)
  keyword <- rep(NA_character_, length(files))
  per_file <- lapply(seq_along(files), function(i) {
    fr <- flowCore::read.FCS(files[i], transformation = FALSE, truncate_max_range = FALSE)
    if (compensate) {
      spill <- Filter(Negate(is.null), flowCore::spillover(fr))
      if (!length(spill))
        stop(basename(files[i]), " has no spillover matrix; pass compensate = FALSE if it is already compensated",
             call. = FALSE)
      keyword[i] <<- names(spill)[1]
      fr <- flowCore::compensate(fr, spill[[1]])
    }
    ch <- fcs_channels(fr, markers)
    ev <- flowCore::exprs(fr)[, ch, drop = FALSE]
    colnames(ev) <- names(ch)
    data.frame(sample_id = sample_id[i], ev, check.names = FALSE)
  })
  out <- do.call(rbind, per_file)
  rownames(out) <- NULL
  attr(out, "compensation") <- data.frame(sample_id = sample_id, file = basename(files), keyword = keyword)
  out
}

#' Read a GatingSet population into a cells table
#'
#' Returns the events of one gated population (the parent, e.g. CD4) of every
#' sample in a flowWorkspace `GatingSet`, in the long `cells` format that
#' [cj_gate()], [cj_from_cells()] and [cj_quickstart()] take. This is the route
#' for data gated in FlowJo (imported with CytoML) or with an openCyto template,
#' when the cytokine gates are to be set by [cj_gate()] rather than read from
#' the GatingSet ([cj_from_gatingset()] reads existing cytokine gates).
#'
#' Values are as stored in the GatingSet, which is usually already compensated
#' and transformed; then use `transform = "none"` in [cj_gate()]. With
#' `inverse_transform = TRUE` the GatingSet's transforms are undone first.
#'
#' @param gs A `GatingSet`.
#' @param node Path or unique name of the population, e.g. `"CD4"`.
#' @param markers As in [cj_read_fcs()]: channel names or exact descriptions,
#'   names give the output columns.
#' @param inverse_transform Undo the GatingSet's transforms
#'   (`flowWorkspace::gs_pop_get_data(inverse.transform = TRUE)`).
#' @return A data frame with `sample_id` (the GatingSet sample names) and one
#'   numeric column per marker.
#' @export
#' @examples
#' \dontrun{
#' gs <- flowWorkspace::load_gs(system.file("extdata", "gs_manual", package = "flowWorkspaceData"))
#' cells <- cj_cells_from_gatingset(gs, "CD4", c(CD38 = "CD38 APC", HLADR = "HLA-DR V500"))
#' }
cj_cells_from_gatingset <- function(gs, node, markers, inverse_transform = FALSE) {
  if (!requireNamespace("flowWorkspace", quietly = TRUE))
    stop("cj_cells_from_gatingset() needs the flowWorkspace package", call. = FALSE)
  markers <- marker_map(markers)
  cs <- flowWorkspace::gs_pop_get_data(gs, node, inverse.transform = inverse_transform)
  ids <- flowWorkspace::sampleNames(gs)
  per_sample <- lapply(ids, function(id) {
    cf <- cs[[id]]
    ch <- fcs_channels(cf, markers)
    ev <- flowCore::exprs(cf)[, ch, drop = FALSE]
    colnames(ev) <- names(ch)
    data.frame(sample_id = id, ev, check.names = FALSE)
  })
  out <- do.call(rbind, per_sample)
  rownames(out) <- NULL
  out
}
