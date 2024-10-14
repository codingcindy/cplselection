#' This function simulates a data set of sample selection with use-specified parameters.
#' 
#' @param outcome_dist Marginal distribution of outcome variable. 
#' @param outcome_par A list of marginal parameters for the outcome equation. Must include "beta".
#' @param select_dist Marginal distribution of selection variable.
#' @param select_par A list of marginal parameters for the selection equation. Must include "beta".
#' @param theta The dependence parameter in the Gaussian copula function.
#' @param x_dist A list of explanatory variable distributional specification lists. 
#'        Options are "Normal", "Uniform", "Poisson", "Negative Binomial" and "Normal Mixture".
#' @param x_par A list of explanatory variable distributional parameters.
#' @param nobs Number of observations.
#' @param selection Whether sample selection happens.
#'
#' @importFrom stats rnorm runif rpois rnbinom pnorm qnorm qpois qnbinom qbinom
#' @export
#' 
simdata <- function(outcome_dist, outcome_par, select_dist, select_par, 
                    theta, nobs, vis, x_dist, selection=TRUE) {
  ##==== get observation counts ====
  ntotal <- round(nobs/vis)
  ##==== generate X variables ====
  ## number of vars
  ko <- length(outcome_par$beta)-1
  ks <- length(select_par$beta)-1
  k  <- max(ko,ks)
  X  <- matrix(NA,nrow=ntotal,ncol=k)
  for (i in (1:k)) {
    dist <- x_dist[[i]]$dist
    if (dist=="Normal") {
      X[,i] <- rnorm(ntotal, mean=x_dist[[i]]$mean, sd=x_dist[[i]]$sd)
    } else if (dist=="Uniform") {
      X[,i] <- runif(ntotal, min=x_dist[[i]]$min, max=x_dist[[i]]$max)
    } else if (dist=="Poisson") {
      X[,i] <- rpois(ntotal, lambda=x_dist[[i]]$lambda)
    } else if (dist=="Negative Binomial") {
      X[,i] <- rnbinom(ntotal, mu=x_dist[[i]]$mu, size=x_dist[[i]]$size)
    } else if (dist=="Normal Mixture") {
      U <- runif(ntotal)
      for (w in (1:length(x_dist[[i]]$weight))) {
        idx <- U<=sum(x_dist[[i]]$weight[1:w])
        idx[is.na(idx)] <- FALSE
        X[idx,i] <- rnorm(n=sum(idx,na.rm=TRUE), mean=x_dist[[i]]$mean[w], sd=x_dist[[i]]$sd[w])
        U[idx] <- NA
      }
    }
  }
  ##==== generate copula variables ====
  ## correlation matrix
  R <- matrix(data=c(1,theta,theta,1), nrow=2)
  ## Gaussian copula variables
  Z <- matrix(data=rnorm(n=ntotal*2), nrow=ntotal)
  C <- t(t(chol(R))%*%t(Z))
  U <- pnorm(C)
  ## ==== generate Y variables ====
  ## ---- generate YS and betaS0 ----
  u <- U[,1]
  xb <- X[,1:ks]%*%select_par$beta[-1]
  y_dist <- select_dist
  b0 <- 0
  ub <- 3
  d0 <- (ub-b0)/2
  if (y_dist=="Probit") {
    y <- qbinom(p=u, prob=pnorm(xb), size=1)
    while (sum(y)!=nobs) {
      if (sum(y)<nobs) { 
        # increase b0
        b0 <- b0+d0
      } else if (sum(y)>nobs) {
        # decrease b0
        b0 <- b0-d0
      } else {}
      d0 <- d0/2
      y <- qbinom(p=u, prob=pnorm(xb+b0), size=1)
    }
  } else if (y_dist=="Logit") {
    y <- qbinom(p=u, prob=1/(1+exp(-xb)), size=1)
    while (sum(y)!=nobs) {
      if (sum(y)<nobs) { 
        # increase b0
        b0 <- b0+d0
      } else if (sum(y)>nobs) {
        # decrease b0
        b0 <- b0-d0
      } else {}
      d0 <- d0/2
      y <- qbinom(p=u, prob=1/(1+exp(-xb-b0)), size=1)
    }
  } else if (y_dist=="Cloglog") {
    y <- qbinom(p=u, prob=1-exp(-exp(xb)), size=1)
    while (sum(y)!=nobs) {
      if (sum(y)<nobs) { 
        # increase b0
        b0 <- b0+d0
      } else if (sum(y)>nobs) {
        # decrease b0
        b0 <- b0-d0
      } else {}
      d0 <- d0/2
      y <- qbinom(p=u, prob=1-exp(-exp(xb+b0)), size=1)
    }
  } else {
    print("Error: Specified selection distribution not supported.")
  }
  ## update intercept & YS
  select_par$beta[1] <- b0
  YS <- y
  ## ---- generate YO ----
  u <- U[,2]
  xb <- cbind(1, X[,1:ko])%*%outcome_par$beta
  y_dist <- outcome_dist
  par <- outcome_par
  if (y_dist=="Probit") {
    y <- qbinom(p=u, prob=pnorm(xb), size=1)
  } else if (y_dist=="Logit") {
    y <- qbinom(p=u, prob=1/(1+exp(-xb)), size=1)
  } else if (y_dist=="Cloglog") {
    y <- qbinom(p=u, prob=1-exp(-exp(xb)), size=1)
  } else if (y_dist=="Normal") {
    y <- qnorm(p=u, mean=xb, sd=par$sigma)
  } else if (y_dist=="Poisson") {
    y <- qpois(p=u, lambda=exp(xb))
  } else if (y_dist=="Negative Binomial") {
    y <- qnbinom(p=u, mu=exp(xb), size=par$r)
  } else if (y_dist=="Exponential") {
    # y <- qbinom(p=u, prob=pnorm(xb), size=1)
  } else {
    print("Error: Specified outcome distribution not supported.")
  }
  YO <- y
  ##==== generate Data Frame ====
  df <- cbind(data.frame(YO=YO, YS=YS), data.frame(X))
  if (selection) {
    ## if selection happends, some observations are missing
    df$YOstar <- df$YO
    df$YO[df$YS==0] <- NA
  }
  res <- list(data=df, intercept=select_par$beta[1])
  
  return(res)
}