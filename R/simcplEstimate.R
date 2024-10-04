#' Simulate sample data and estimate Bayesian copula sample selection model
#'
#' @param seed 
#' @param filedir 
#' @param filename 
#' @param outcome_dist 
#' @param select_dist 
#' @param outcome_par 
#' @param select_par 
#' @param theta 
#' @param x_dist 
#' @param nobs 
#' @param outcome_formula 
#' @param select_formula 
#' @param loop 
#' @param burnin 
#' @param stepsize 
#' @param stepadj 
#' @param sliceadj 
#'
#' @return NULL
#' @export 
#'
#' @examples
simcplEstimate <- function(
    seed, filedir, filename,
    outcome_dist, select_dist, outcome_par, select_par, theta, x_dist, nobs, 
    outcome_formula, select_formula, loop, burnin, stepsize, stepadj, sliceadj
    ) {
  ## simulate data set
  set.seed(seed)
  df <- simdata(outcome_dist=outcome_dist, outcome_par=outcome_par,
                select_dist=select_dist, select_par=select_par,
                theta=theta, x_dist=x_dist, nobs=nobs, selection=TRUE)
  ## get true values
  if (outcome_dist=="Normal") {
    morepar <- outcome_par$sigma
  } else if (outcome_dist=="Negative Binomial") {
    morepar <- outcome_par$r
  } else {
    morepar <- NULL
  }
  trueval <- c(theta=theta, betaO=outcome_par$beta[1:3], 
               betaS=select_par$beta, morepar=morepar)
  ## copula sample selection model
  set.seed(seed)
  iterates <- cplselectionMCMC(outcome_formula, select_formula,
                               outcome_dist, select_dist, data=df,
                               loop, stepsize, stepadj, sliceadj)
  res <- cplselectionInfer(iterates, burnin)
  ## data to save
  reslist <- list()
  reslist$outcome_dist <- outcome_dist
  reslist$select_dist <- select_dist
  reslist$x_dist <- x_dist
  reslist$outcome_par <- outcome_par
  reslist$select_par <- select_par
  reslist$theta <- theta
  reslist$morepar <- morepar
  iterates$ZS <- NULL
  iterates$MH_step <- NULL
  iterates$acceptance_prob <- NULL
  reslist$iterates <- iterates
  reslist$result <- res
  ## save results
  fname <- paste0(filedir,"/",filename,"_seed",seed,".RData")
  save(reslist, file=fname)
  ## ensure the object size
  print(format(object.size(reslist), units="MB"))
  return(NULL)
}