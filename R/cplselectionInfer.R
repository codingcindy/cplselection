#' Conduct Bayesian inference with given MCMC iterates
#'
#' @param iterates A list of MCMC iterates matrix returned by `cplselectionMCMC` funciton
#' @param burnin Burn-in iterations to discard.
#'
#' @importFrom stats sd quantile
#' @return A data.frame object containing estimation results.
#' @export
#'
#' @examples NULL
cplselectionInfer <- function(iterates, burnin) {
  loop <- dim(iterates$ZS)[2]
  range <- burnin:loop
  iterates$ZS <- NULL
  iterates$MH_step <- NULL
  iterates$MH_acceptance <- NULL
  iterates$MH_acceptance_prob <- NULL
  res <- data.frame()
  for (est in names(iterates)) {
    if (is.null(dim(iterates[[est]]))) {
      iters <- iterates[[est]][range]
      rescoef <- matrix(c(sd(iters),quantile(iters,probs=c(0.5,0.025,0.975))),nrow=1)
      colnames(rescoef) <- c("sd","median","lb","ub")
      rownames(rescoef) <- est
    } else {
      iters <- iterates[[est]][,range]
      rescoef <- t(rbind(
        apply(iters, MARGIN=1, FUN=sd),
        apply(iters, MARGIN=1, FUN=quantile, probs=c(0.5,0.025,0.975))
        ))
      colnames(rescoef) <- c("sd","median","lb","ub")
      rownames(rescoef) <- paste0(est,0:(nrow(iters)-1))
    }
    res <- rbind(res,rescoef)
  }
  res <- res[,c("median","sd","lb","ub")]
  return(res)
}