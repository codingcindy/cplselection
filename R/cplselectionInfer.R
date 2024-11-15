#' Infer estimates from MCMC iterates.
#' 
#' This function makes Bayesian inference with given MCMC iterates returned by `cplselectionMCMC` function.
#'
#' @param iterates A list of MCMC iterates returned by `cplselectionMCMC` function.
#' @param burnin Burn-in iterations to discard.
#' @param trueval True parameters, if available. Default set to NULL.
#'
#' @importFrom stats sd quantile
#' @return A data.frame object containing estimation results.
#'
#' @examples 
#' # Example data : Poisson outcome variable YO, logit selection process YS
#' nobs <- 300
#' dependence <- 0.5
#' betaO0 <- 1
#' betaO1 <- 1
#' betaS0 <- 0
#' betaS1 <- 1
#' betaS2 <- 1
#' uo <- runif(nobs)
#' us <- pnorm(dependence*qnorm(uo) + dependence*rnorm(nobs))
#' data <- data.frame(X1 = rnorm(nobs), X2 = rnorm(nobs))
#' data$YO <- qpois(p=uo, lambda=exp(betaO0 + betaO1*data$X1))
#' data$YS <- qbinom(p=us, prob=1/(1+exp(-(betaS0 + betaS1*data$X1 + betaS2*data$X2))), size=1)
#' data$YO[data$YS==0] <- NA
#' 
#' # Perform MCMC draws
#' draws <- cplselectionMCMC(outcome_formula = YO~X1, select_formula = YS~X1+X2,
#'                           outcome_dist = "Poisson", select_dist = "Logit", 
#'                           data = data, loop = 5000)
#'                           
#' # Conduct Inference
#' # If the true values are not known, returns estimates:
#' result.1 <- cplselectionInfer(iterates=draws, burnin=loop/2) 
#' print(round(result.1, 2))
#' 
#' # If the true values are known, returns true values together with estimates:
#' trueval <- c(dependence, betaO0, betaO1, betaS0, betaS1, betaS2)
#' result.2 <- cplselectionInfer(iterates=draws, burnin=loop/2, trueval=trueval)
#' print(round(result.2, 2))
#' 
#' @export
cplselectionInfer <- function(iterates, burnin, trueval=NULL) {
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
  if (!is.null(trueval)) {
    res <- cbind(trueval, res)
  }
  return(res)
}