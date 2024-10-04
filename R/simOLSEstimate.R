#' Perform OLS estimation at scale
#'
#' @param seed 
#' @param outcome_dist 
#' @param select_dist 
#' @param outcome_par 
#' @param select_par 
#' @param theta 
#' @param x_dist 
#' @param nobs 
#' @param outcome_formula 
#' @param select_formula 
#'
#' @importFrom utils head 
#' @importFrom stats lm glm binomial
#' @importFrom MASS glm.nb
#' 
#' @return A data.frame object containing OLS estimation results.
#' @export
#'
#' @examples NULL
simOLSEstimate <- function(
    seed, outcome_dist, select_dist, outcome_par, select_par, theta, x_dist, nobs, 
    outcome_formula, select_formula) {
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
  trueval <- c(betaO=head(outcome_par$beta,-1), morepar=morepar)
  cnames <- paste0("betaO", 0:(length(outcome_par$beta)-2))
  ## naive OLS model
  if (outcome_dist=="Normal") {
    tmp <- summary(lm(formula=outcome_formula, data=df))
    res <- rbind(tmp$coefficients, tmp$sigma)
    cnames <- c(cnames,"sigmaO")
  } else if (outcome_dist=="Probit") {
    res <- summary(glm(outcome_formula, family=binomial(link="probit"), data=df))$coefficients
  } else if (outcome_dist=="Logit") {
    res <- summary(glm(outcome_formula, family=binomial(link="logit"), data=df))$coefficients
  } else if (outcome_dist=="Poisson") {
    res <- summary(glm(outcome_formula, family="poisson", data=df))$coefficients
  } else if (outcome_dist=="Negative Binomial") {
    tmp <- summary(MASS::glm.nb(outcome_formula, data=df))
    res <- rbind(tmp$coefficients, c(tmp$theta, tmp$SE.theta, NA, NA))
    cnames <- c(cnames,"rO")
  } else {
    stop("Outcome distribution not supported.")
  }
  ## data to save
  resdf <- data.frame(
    outcome_dist=outcome_dist, 
    select_dist=select_dist,
    x_dist=x_dist[[1]]$dist,
    dependence=theta,
    selectivity=ifelse((x_dist[[1]]$dist=="Poisson"),
                       select_par$beta[1]+1, select_par$beta[1]),
    pctobs=mean(df$YS),
    seed=seed,
    parameter=cnames,
    trueval=trueval,
    estimate=res[,1],
    se=res[,2],
    pval=res[,4]
  )
  
  return(resdf)
}