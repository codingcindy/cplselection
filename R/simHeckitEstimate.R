#' Simulate data and estimate with Heckit method
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
#' @return
#' @export
#'
#' @examples
simHeckitEstimate <- function(
    seed, outcome_dist, select_dist, outcome_par, select_par, theta, x_dist, nobs, 
    outcome_formula, select_formula) {
  library(MASS)
  library(stats)
  library(sampleSelection)
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
  trueval <- c(betaS=select_par$beta, betaO=head(outcome_par$beta,-1), 
               IMR=NA, morepar=morepar, rho=theta)
  cnames <- c(paste0("betaS", 0:(length(select_par$beta)-1)), 
              paste0("betaO", 0:(length(outcome_par$beta)-2)), 
              "IMR")
  ## heckit sample selection model
  if (outcome_dist=="Normal") {
    res <- summary(heckit(selection=select_formula, outcome=outcome_formula, 
                          data=df, method="2step"))$estimate
    cnames <- c(cnames,"sigmaO","rho")
  } else {
    ## first stage
    step1 <- glm(select_formula, family=binomial(link="probit"), data=df)
    res1 <- summary(step1)$coefficients
    xb1 <- predict(step1)
    IMR <- dnorm(xb1)/pnorm(xb1)
    df$IMR <- IMR
    ## second stage
    vec <- all.vars(outcome_formula)
    outcome_eq <- reformulate(c(vec[2:length(vec)], "IMR"), vec[1])
    if (outcome_dist=="Probit") {
      step2 <- glm(outcome_eq, family=binomial(link="probit"), data=df)
      cnames <- c(cnames,"rho")
    } else if (outcome_dist=="Logit") {
      step2 <- glm(outcome_eq, family=binomial(link="logit"), data=df)
      cnames <- c(cnames,"rho")
    } else if (outcome_dist=="Poisson") {
      step2 <- glm(outcome_eq, family="poisson", data=df)
      cnames <- c(cnames,"rho")
    } else if (outcome_dist=="Negative Binomial") {
      step2 <- MASS::glm.nb(outcome_eq, data=df)
      # overdispersion
      morepar <- c(step2$theta,NA,NA,NA)
      cnames <- c(cnames,"thetaO","rho")
    } else {
      stop("Outcome distribution not supported.")
    }
    res2 <- summary(step2)$coefficients
    res <- rbind(res1, res2, morepar=morepar, rho=NA)
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