# Hand-built data: 40 donors, visits pre and post, one stimulation. Positive cells sit in
# the IFNG-only pattern. The first `n_resp` donors respond at post; the rest do not.
mimosa2_data <- function(n_donor = 40, n_resp = 20, N = 20000L, seed = 3) {
  set.seed(seed)
  bits <- cj_bits(c("IFNG", "TNF"))
  lab <- pattern_labels(bits)
  ifng_only <- which(bits[, "IFNG"] == 1L & bits[, "TNF"] == 0L)
  neg <- which(rowSums(bits) == 0L)
  donors <- sprintf("D%02d", seq_len(n_donor))
  samples <- expand.grid(condition = c("control", "CMV"), visit = c("pre", "post"), donor_id = donors,
                         stringsAsFactors = FALSE)
  samples$visit <- factor(samples$visit, levels = c("pre", "post"))
  samples$parent_population <- "CD4"
  samples$sample_id <- paste(samples$donor_id, samples$visit, samples$condition, sep = "_")
  samples$N <- N
  resp <- samples$donor_id %in% donors[seq_len(n_resp)]
  p <- ifelse(samples$condition == "control", 0.001,
              ifelse(samples$visit == "post" & resp, 0.006, 0.002))
  counts <- matrix(0L, nrow(samples), length(lab), dimnames = list(samples$sample_id, lab))
  counts[, ifng_only] <- as.integer(rbinom(nrow(samples), N, p))
  counts[, neg] <- N - counts[, ifng_only]
  stim <- samples[samples$condition == "CMV", ]
  contrasts <- data.frame(contrast_id = paste(stim$donor_id, stim$visit, "CD4_CMV", sep = "_"),
                          donor_id = stim$donor_id, visit = stim$visit, parent_population = "CD4",
                          stimulation = "CMV", stim_sample_id = stim$sample_id,
                          control_sample_id = sub("_CMV$", "_control", stim$sample_id))
  samples <- samples[c("sample_id", "donor_id", "visit", "parent_population", "condition", "N")]
  list(x = cytoJointR:::new_cj_data(counts, samples, contrasts, bits), responders = donors[seq_len(n_resp)])
}
