# Lightweight independent simulation of paired cell-level data with known
# multinomial pattern probabilities and threshold-compatible fluorescence.
# Stress generators with WRONG emission models must be added for benchmarks.
#' @rdname cytoJointR-reference
#' @export
cj_simulate_cells <- function(donors=12L, cells_per_sample=1500L,
                              markers=c("IFNG", "TNF", "IL2"),
                              regime=c("antigen", "polyclonal"),
                              missing_controls=0.2, seed=1L,
                              threshold=0) {
  regime <- match.arg(regime)
  stopifnot(donors>=2, cells_per_sample>=1,
            length(markers)>=2, length(markers)<=5,
            missing_controls>=0, missing_controls<1)
  set.seed(seed)
  bits <- cj_bits(markers)
  J <- length(markers); K <- nrow(bits)
  baseline_intercept <- if (regime=="antigen") -5.8 else -3.5
  stimulus_effect <- if (regime=="antigen") 1.8 else 4.2
  d <- rowSums(bits)
  pair <- (d*(d-1))/2
  control_missing <- sample.int(donors,size=floor(donors*missing_controls))
  raw_cells <- vector("list", donors*2L)
  truth <- vector("list",donors)
  n_out <- 0L
  for (i in seq_len(donors)) {
    grp <- as.integer(i > floor(donors/2))
    donor_shift <- rnorm(1,0,0.35)
    p_condition <- vector("list",2L)
    for (stim in 0:1) {
      eta <- d*(baseline_intercept + donor_shift +
                     stim*(stimulus_effect + 0.35*grp)) +
             pair*(0.15 + stim*0.25)
      eta <- eta-max(eta)
      prob <- exp(eta)/sum(exp(eta))
      p_condition[[stim+1L]] <- prob
      if (stim==0L && i %in% control_missing) next
      codes <- sample.int(K,cells_per_sample,replace=TRUE,prob=prob)
      Z <- bits[codes,,drop=FALSE]
      X <- matrix(0,cells_per_sample,J)
      for (j in seq_len(J)) {
        mu <- ifelse(Z[,j]==1L,1.6,-1.6) + stim*0.2
        # Sample EXACTLY from the Gaussian conditional on the sign-defined gate.
        pcut <- pnorm(threshold,mean=mu,sd=0.8)
        u <- runif(cells_per_sample)
        cdf <- ifelse(Z[,j]==1L,pcut+(1-pcut)*u,pcut*u)
        cdf <- pmax(1e-12,pmin(1-1e-12,cdf))
        X[,j] <- qnorm(cdf,mean=mu,sd=0.8)
      }
      colnames(X) <- markers
      n_out <- n_out+1L
      raw_cells[[n_out]] <- data.frame(
        sample_id=sprintf("D%03d_%s",i,if (stim) "stim" else "unstim"),
        donor_id=sprintf("D%03d",i),
        condition=if (stim) "stim" else "unstim",
        cell_subset="CD4_T", stimulus="test_stimulation",
        group=grp, X, check.names=FALSE
      )
    }
    score <- cj_scores_from_probs(p_condition[[1]],p_condition[[2]],bits)
    truth[[i]] <- data.frame(donor_id=sprintf("D%03d",i),group=grp,
                            has_control=!(i %in% control_missing),
                            t(score), check.names=FALSE)
  }
  list(cells=do.call(rbind,raw_cells[seq_len(n_out)]),
       truth=do.call(rbind,truth),bits=bits)
}
