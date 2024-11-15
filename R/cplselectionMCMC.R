#' MCMC sampling of the Bayesian copula sample selection model
#' 
#' This function performs MCMC sampling of the Bayesian copula sample selection model using Gibbs within Metropolis-Hastings algorithm.
#' 
#' @param outcome_formula The outcome equation formula, e.g. YO ~ X1.
#' @param select_formula The selection equation formula, e.g. YS ~ X1 + X2.
#' @param outcome_dist The outcome variable distribution, supporting binomial with probit link ("Probit"), binomial with logit link ("Logit"), binomial with complementary log-log link ("Cloglog"), Poisson ("Poisson"), negative binomial ("Negative Binomial"), exponential ("Exponential") and normal ("Normal") distributions. 
#' @param select_dist The selection variable distribution, supporting binomial distribution with probit link ("Probit"), binomial with logit link ("Logit"), and binomial with complementary log-log link ("Cloglog").
#' @param data The estimation sample data set, can be a data.frame object or a named matrix with column as variables.
#' @param loop The maximum number of iterations for MCMC draws.
#' @param stepsize The step size in the Metropolis Hastings algorithm, refers to the magnitude of proposed moves in the parameter space.
#' @param stepadj The adjustment proportion of the MH step size based on the acceptance rates in the current slice.
#' @param sliceadj The number of draws in a step adjustment slice.
#' 
#' @importFrom progress progress_bar
#' @importFrom truncnorm rtruncnorm
#' @importFrom MASS mvrnorm
#' @importFrom stats runif rnorm pnorm qnorm ppois pnbinom terms
#' 
#' @returns A list of iteration results, including MCMC parameter draws, acceptance probabilities and step sizes per slice.
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
#' @export
cplselectionMCMC <- function(outcome_formula, select_formula, 
                             outcome_dist, select_dist, data, 
                             loop, stepsize=0.1, stepadj=0.05, 
                             sliceadj=500) {
  ##==== Auxiliary Functions ====
  logPosteriorRho <- function(rho,zs,zo) {
    if (abs(rho)>=1) {
      log_posterior <- -Inf 
    } else {
      log_posterior <- (-length(zo)/2)*log(1-rho^2)+
        (-1/2)*sum(zo^2-2*rho*zs*zo+zs^2)/(1-rho^2)
    }
    return(log_posterior)
  }
  drawBeta <- function(y,x,sigma,zcond,rho,V0inv,V0invu0) {
    Omegainvdiag <- 1/(1-rho^2)
    Omegainvoffdiag <- rho/(rho^2-1)
    Vn <- solve(V0inv+Omegainvdiag*t(x)%*%as.matrix(x)/(sigma^2))
    Vn <- (Vn+t(Vn))/2
    un <- Vn%*%(V0invu0+Omegainvdiag*t(x)%*%as.matrix(y)/(sigma^2)+Omegainvoffdiag*t(x)%*%zcond/sigma)
    draw <- MASS::mvrnorm(1,un,Vn)
    return(draw)
  }
  logPosteriorBeta <- function(beta,y,x,zcond,rho,u0,V0inv,y_dist,morepar=NULL) {
    # log prior density
    log_prior <- (-1/2)*(beta-u0)%*%V0inv%*%(beta-u0)
    # log likelihood
    if (y_dist=="Probit") {
      cutoff <- -x%*%beta
      prob <- pnorm((cutoff-rho*zcond)/sqrt(1-rho^2))
      prob[is.na(zcond)] <- pnorm(cutoff[is.na(zcond)])
      prob[prob==0] <- .Machine$double.xmin
      prob[prob==1] <- 1-1e+300*.Machine$double.xmin
      log_lik <- sum((1-y)*log(prob)+y*log(1-prob))
    } else if (y_dist=="Logit") {
      cutoff <- qnorm(1-1/(1+exp(-x%*%beta)))
      prob <- pnorm((cutoff-rho*zcond)/sqrt(1-rho^2))
      prob[is.na(zcond)] <- pnorm(cutoff[is.na(zcond)])
      prob[prob==0] <- .Machine$double.xmin
      prob[prob==1] <- 1-1e+300*.Machine$double.xmin
      log_lik <- sum((1-y)*log(prob)+y*log(1-prob))
    } else if (y_dist=="Cloglog") {
      cutoff <- qnorm(exp(-exp(x%*%beta)))
      prob <- pnorm((cutoff-rho*zcond)/sqrt(1-rho^2))
      prob[is.na(zcond)] <- pnorm(cutoff[is.na(zcond)])
      prob[prob==0] <- .Machine$double.xmin
      prob[prob==1] <- 1-1e+300*.Machine$double.xmin
      log_lik <- sum((1-y)*log(prob)+y*log(1-prob))
    } else if (y_dist=="Poisson") {
      lambda <- exp(x%*%beta)
      TU <- qnorm(ppois(y,lambda,log.p=TRUE),log.p=TRUE)
      TL <- qnorm(ppois(y-1,lambda,log.p=TRUE),log.p=TRUE)
      lik <- rep(0,length(zcond))
      lik[!is.na(zcond)] <- pnorm((TU[!is.na(zcond)]-rho*zcond[!is.na(zcond)])/sqrt(1-rho^2))-
        pnorm((TL[!is.na(zcond)]-rho*zcond[!is.na(zcond)])/sqrt(1-rho^2))
      lik[is.na(zcond)] <- pnorm(TU[is.na(zcond)])-pnorm(TL[is.na(zcond)])
      lik[lik==0] <- .Machine$double.xmin
      lik[lik==1] <- 1-1e+300*.Machine$double.xmin
      log_lik <- sum(log(lik))
    } else if (y_dist=="Negative Binomial") {
      if (is.null(morepar)) {
        stop("Provide a probability parameter for Negative Binomial distribution.")
      }
      mu <- exp(x%*%beta)
      TU <- qnorm(pnbinom(q=y, mu=mu, size=morepar, log.p=TRUE), log.p=TRUE)
      TL <- qnorm(pnbinom(q=y-1, mu=mu, size=morepar, log.p=TRUE), log.p=TRUE)
      log_lik <- sum(log(
        pnorm((TU-rho*zcond)/sqrt(1-rho^2)) - pnorm((TL-rho*zcond)/sqrt(1-rho^2)))
        )
    } else if (y_dist=="Exponential") {
      xb <- x%*%beta
      rate <- exp(xb)
      z <- qnorm(1-exp(-rate*y))
      log_lik <- sum(xb-rate*y-0.5*(rho^2)*(z^2)/(1-rho^2)+zcond*z*rho/(1-rho^2))
    } else {
      stop("Specified marginal distribution not supported.")
    }
    log_posterior <- log_prior+log_lik
    return(log_posterior)
  }
  logPosteriorSigma <- function(sigma,y,x,beta,zcond,rho,a0,b0) {
    nobs <- length(y)
    log_prior <- -(a0+2+nobs)*log(sigma)
    OmegainvDiag <- 1/(1-rho^2)
    OmegainvOffdiag <- rho/(rho^2-1)
    z <- (y-x%*%beta)/sigma
    log_lik <- (-1/2)*(b0/sigma^2+OmegainvDiag*sum(z^2)+2*OmegainvOffdiag*sum(zcond*z))
    log_posterior <- log_prior+log_lik
    return(log_posterior)
  }
  logPosteriorR <- function(r,y,x,beta,zcond,rho) {
    log_prior <- exp(-r)
    exb <- exp(x%*%beta)
    TU <- qnorm(pnbinom(y,mu=exb,size=r,log.p=TRUE),log.p=TRUE)
    TL <- qnorm(pnbinom(y-1,mu=exb,size=r,log.p=TRUE),log.p=TRUE)
    lik <- pnorm((TU-rho*zcond)/sqrt(1-rho^2))-pnorm((TL-rho*zcond)/sqrt(1-rho^2))
    lik[lik==0] <- .Machine$double.xmin
    lik[lik==1] <- 1-1e+300*.Machine$double.xmin
    log_lik <- sum(log(lik))
    log_posterior <- log_prior+log_lik
    return(log_posterior)
  }
  computezCont <- function(y_dist,y,x,beta,morepar=NULL) {
    if (y_dist=="Normal") {
      if (is.null(morepar)) {
        stop("Provide variance parameter for Normal distribution.")
      }
      res <- (y-x%*%beta)/morepar
    } else if (y_dist=="Exponential") {
      res <- qnorm(pexp(q=y, rate=exp(x%*%beta), log.p=TRUE), log.p=TRUE)
    } else {
      ## return error message and stop execution
      stop("Specified distribution not supported.")
    }
  }
  drawzDisc <- function(y_dist,y,x,beta,zcond,rho,morepar=NULL) {
    if (y_dist=="Probit") {
      cutoff <- -x%*%beta
      TU <- cutoff
      TL <- cutoff
      TU[y==1] <- Inf
      TL[y==0] <- -Inf
    } else if (y_dist=="Logit") {
      cutoff <- qnorm(1-1/(1+exp(-x%*%beta)))
      TU <- cutoff
      TL <- cutoff
      TU[y==1] <- Inf
      TL[y==0] <- -Inf
    } else if (y_dist=="Cloglog") {
      cutoff <- qnorm(exp(-exp(x%*%beta)))
      TU <- cutoff
      TL <- cutoff
      TU[y==1] <- Inf
      TL[y==0] <- -Inf
    } else if (y_dist=="Poisson") {
      lambda <- exp(x%*%beta)
      TU <- qnorm(ppois(y,lambda,log.p=TRUE),log.p=TRUE)
      TL <- qnorm(ppois(y-1,lambda,log.p=TRUE),log.p=TRUE)
    } else if (y_dist=="Negative Binomial") {
      if (is.null(morepar)) {
        stop("Provide a probability parameter for Negative Binomial distribution.")
      }
      mu <- exp(x%*%beta)
      TU <- qnorm(pnbinom(y,mu=mu,size=morepar,log.p=TRUE),log.p=TRUE)
      TL <- qnorm(pnbinom(y-1,mu=mu,size=morepar,log.p=TRUE),log.p=TRUE)
    } else {
      stop("Specified distribution not supported.")
    }
    TL[TL==Inf] <- .Machine$double.xmax/10
    mean_cond <- rho*zcond
    mean_cond[is.na(mean_cond)] <- 0
    sd_cond <- rep(sqrt(1-rho^2), length(zcond))
    sd_cond[is.na(mean_cond)] <- 1
    draws <- truncnorm::rtruncnorm(n=1,a=TL,b=TU,mean=mean_cond,sd=sd_cond)
    return(draws)
  }
  
  ##==== Sampling Setups ====
  ##---- Dimensions and design matrices ----
  XO <- as.matrix(data[,all.vars(outcome_formula[[3]])])
  XS <- as.matrix(data[,all.vars(select_formula[[3]])])
  if (attr(terms(outcome_formula),"intercept")==1) {
    XO <- cbind(1, XO)
  }
  if (attr(terms(select_formula),"intercept")==1) {
    XS <- cbind(1, XS)
  }
  YO <- as.matrix(eval(outcome_formula[[2]], envir=data))
  YS <- as.matrix(eval(select_formula[[2]], envir=data))
  nobs <- nrow(data)
  KO <- ncol(XO)
  KS <- ncol(XS)
  ##---- Iteration containers & Initiation ----
  parnames <- c("rho", paste0("betaO",0:(KO-1)), paste0("betaS",0:(KS-1)))
  MH_par <- parnames
  zo <- rep(NA, nobs)
  ZS <- cbind(0, matrix(NA, nobs, loop))
  rho <- cbind(0, matrix(NA, 1, loop))
  betaO <- cbind(0, matrix(NA, KO, loop))
  betaS <- cbind(0, matrix(NA, KS, loop))
  moreparO <- NULL
  moreparS <- NULL
  for (item in (1:2)) {
    y_dist <- c(outcome_dist, select_dist)[item]
    suffix <- c("O", "S")[item]
    if (y_dist=="Normal") {
      assign(parnam <- paste0("sigma", suffix), 
             cbind(1, matrix(NA, 1, loop)))
      parnames <- c(parnames, parnam)
      MH_par <- setdiff(parnames, paste0("beta",suffix,0:(KS-1)))
      assign(paste0("morepar", suffix), parnam)
    } else if (y_dist=="Negative Binomial") {
      assign(parnam <- paste0("r", suffix), 
             cbind(0.5, matrix(NA, 1, loop)))
      parnames <- c(parnames, parnam)
      assign(paste0("morepar", suffix), parnam)
      MH_par <- c(MH_par, parnam)
    } else if (y_dist=="Probit"|y_dist=="Logit"|y_dist=="Cloglog"|
               y_dist=="Poisson"|y_dist=="Exponential") {
      ## no extra parameters required
    } else {
      ## return error message and stop execution
      stop("Specified distribution not supported.")
    }
  }
  if (outcome_dist=="Normal") {
    zo[YS==1] <- computezCont(outcome_dist,y=YO[YS==1],x=XO[YS==1,],beta=betaO[,1],morepar=sigmaO[,1])
  } else if (outcome_dist=="Probit"|outcome_dist=="Logit"|outcome_dist=="Cloglog"|outcome_dist=="Poisson") {
    zo[YS==1] <- drawzDisc(outcome_dist,YO[YS==1],XO[YS==1,],betaO[,1],ZS[YS==1,1],rho[,1])
  } else if (outcome_dist=="Negative Binomial") {
    zo[YS==1] <- drawzDisc(outcome_dist,YO[YS==1],XO[YS==1,],betaO[,1],ZS[YS==1,1],rho[,1],morepar=rO[,1])
  } else if (outcome_dist=="Exponential") {
    zo[YS==1] <- computezCont(outcome_dist,y=YO[YS==1],x=XO[YS==1,],beta=betaO[,1])
  } else {
    stop("Specified distribution not supported.")
  }
  
  ##---- MH acceptance & step settings ----
  MH_n <- length(MH_par)
  MH_step <- matrix(stepsize, nrow=floor(loop/sliceadj), ncol=MH_n)
  MH_acceptance <- data.frame(matrix(0, nrow=loop, ncol=MH_n))
  MH_acceptance_prob <- data.frame(matrix(0, nrow=loop, ncol=MH_n))
  colnames(MH_step) <- MH_par
  colnames(MH_acceptance) <- MH_par
  colnames(MH_acceptance_prob) <- MH_par
  ##==== Priors ====
  ## for betaO
  u0O         <- rep(0,KO)
  V0invO      <- diag(KO)
  V0invu0O    <- V0invO%*%u0O
  ## for betaS
  u0S         <- rep(0,KS)
  V0invS      <- diag(KS)
  V0invu0S    <- V0invS%*%u0S
  ## for sigmaO
  a0          <- 2.5
  b0          <- 1.5
  ## for r2
  # r2 ~ exp(1)
  
  ##==== Iteration Steps ====
  s <- 1;
  pb <- progress::progress_bar$new(
    format="    sampling progress [:bar] :percent in :elapsed | eta: :eta", 
    total=loop, width=80, clear=FALSE)
  for(i in 2:(loop+1)){
    pb$tick()
    ##---- adjust step sizes ----
    if ((i%%sliceadj==1)&(i>sliceadj)&(s<nrow(MH_step))) {
      s <- s+1
      ratio <- colMeans(MH_acceptance[(i-sliceadj-1):(i-1),])
      duminc <- (ratio>0.35)
      dumdec <- (ratio<0.2)
      MH_step[s,] <- MH_step[s-1,]*(1+stepadj*(duminc-dumdec))
    }
    ##---- copula dependence ---- 
    ## update rho: MH with normal proposal ----
    # print("update rho")
    current <- rho[,i-1]
    propose <- truncnorm::rtruncnorm(1,mean=current,sd=MH_step[s,"rho"],a=-1,b=1)
    score_c <- logPosteriorRho(current,ZS[YS==1,i-1],zo[YS==1])
    score_p <- logPosteriorRho(propose,ZS[YS==1,i-1],zo[YS==1])
    if (is.na(score_p-score_c)) {
      logA  <- -Inf
    } else {
      logA  <- min(0, score_p-score_c)
    }
    MH_acceptance_prob[i,"rho"] <- logA
    if (logA > log(runif(1))) {
      rho[,i] <- propose
      MH_acceptance[i,"rho"] <- 1
    } else {
      rho[,i] <- current
    }
    ##---- outcome margin ----
    # print("update betaO")
    ## update betaO: ----
    ##    Gibbs for normal, 
    ##    MH with normal proposal for else
    if (outcome_dist=="Normal") {
      betaO[,i] <- drawBeta(y=YO[YS==1],x=XO[YS==1,],sigma=get(moreparO)[,i-1],
                            zcond=ZS[YS==1,i-1],rho=rho[,i],V0inv=V0invO,V0invu0=V0invu0O)
      # betaO[,i] <- drawBeta(y=YO,x=XO,sigma=get(moreparO)[,i-1],zcond=ZS[,i-1],rho=rho[,i],V0inv=V0invO,V0invu0=V0invu0O)
    } else {
      cnames <- paste0("betaO",0:(KO-1))
      current <- betaO[,i-1]
      propose <- MASS::mvrnorm(1,mu=current,Sigma=diag(MH_step[s,cnames]))
      if (outcome_dist=="Negative Binomial") {
        morepar <- get(moreparO)[,i-1]
      } else {
        morepar <- NULL
      }
      score_c <- logPosteriorBeta(beta=current,y=YO[YS==1],x=XO[YS==1,],zcond=ZS[YS==1,i-1],
                                  rho=rho[,i],u0=u0O,V0inv=V0invO,outcome_dist,morepar)
      score_p <- logPosteriorBeta(beta=propose,y=YO[YS==1],x=XO[YS==1,],zcond=ZS[YS==1,i-1],
                                  rho=rho[,i],u0=u0O,V0inv=V0invO,outcome_dist,morepar)
      if (is.na(score_p-score_c)) {
        logA  <- -Inf
      } else {
        logA  <- min(0, score_p-score_c)
      }
      MH_acceptance_prob[i,cnames] <- logA
      if (logA > log(runif(1))) {
        betaO[,i] <- propose
        MH_acceptance[i,cnames] <- 1
      } else {
        betaO[,i] <- current
      }
    }
    # print("update moreparO")
    ## update other parameters: ----
    ##    MH with log-normal proposal for normal, 
    ##    MH with beta proposal for negative binomial,
    ##    none for else
    if (outcome_dist=="Normal") {
      current <- sigmaO[,i-1]
      propose <- sqrt(exp(rnorm(1,mean=2*log(current),sd=MH_step[s,"sigmaO"])))
      score_c <- logPosteriorSigma(current,YO[YS==1],XO[YS==1,],betaO[,i],ZS[YS==1,i-1],rho[,i],a0,b0)
      score_p <- logPosteriorSigma(propose,YO[YS==1],XO[YS==1,],betaO[,i],ZS[YS==1,i-1],rho[,i],a0,b0)
      if (is.na(score_p-score_c+2*log(propose)-2*log(current))) {
        logA  <- -Inf
      } else {
        logA  <- min(0, score_p-score_c+2*log(propose)-2*log(current))
      }
      MH_acceptance_prob[i,"sigmaO"] <- logA
      if (logA > log(runif(1))) {
        sigmaO[,i] <- propose
        MH_acceptance[i,"sigmaO"] <- 1
      } else {
        sigmaO[,i] <- current
      }
    } else if (outcome_dist=="Negative Binomial") {
      current <- rO[,i-1]
      propose <- exp(rnorm(1,mean=log(current),sd=MH_step[s,"rO"]))
      score_c <- logPosteriorR(r=current,y=YO[YS==1],x=XO[YS==1,],beta=betaO[,i],zcond=ZS[YS==1,i-1],rho=rho[,i])
      score_p <- logPosteriorR(r=propose,y=YO[YS==1],x=XO[YS==1,],beta=betaO[,i],zcond=ZS[YS==1,i-1],rho=rho[,i])
      if (is.na(score_p-score_c+2*log(propose)-2*log(current))) {
        logA  <- -Inf
      } else {
        logA  <- min(0, score_p-score_c+2*log(propose)-2*log(current))
      }
      MH_acceptance_prob[i,"rO"] <- logA
      if (logA > log(runif(1))) {
        rO[,i] <- propose
        MH_acceptance[i,"rO"] <- 1
      } else {
        rO[,i] <- current
      }
    } else {
      # no more parameters to estimate 
    }
    # print("update zo")
    ## update zo: ----
    if (outcome_dist=="Normal") {
      zo[YS==1] <- computezCont(outcome_dist,y=YO[YS==1],x=XO[YS==1,],beta=betaO[,i],morepar=sigmaO[,i])
    } else if (outcome_dist=="Negative Binomial") { 
      zo[YS==1] <- drawzDisc(outcome_dist,y=YO[YS==1],x=XO[YS==1,],beta=betaO[,i],zcond=ZS[YS==1,i-1],rho=rho[,i],morepar=rO[,i])
    } else if (outcome_dist=="Probit"|outcome_dist=="Logit"| outcome_dist=="Cloglog"|outcome_dist=="Poisson") {
      zo[YS==1] <- drawzDisc(outcome_dist,y=YO[YS==1],x=XO[YS==1,],beta=betaO[,i],zcond=ZS[YS==1,i-1],rho=rho[,i])
    } else if (outcome_dist=="Exponential") {
      zo[YS==1] <- computezCont(outcome_dist,y=YO[YS==1],x=XO[YS==1,],beta=betaO[,i])
    } else {
      stop("Specified distribution not supported.")
    }
    
    ##---- selection margin ----
    # print("update betaS")
    ## update betaS: MH with normal proposal ----
    cnames <- paste0("betaS",0:(KS-1))
    current <- betaS[,i-1]
    propose <- MASS::mvrnorm(1,mu=current,Sigma=diag(MH_step[s,cnames]))
    score_c <- logPosteriorBeta(beta=current,y=YS,x=XS,zcond=zo,rho=rho[,i],u0S,V0invS,y_dist=select_dist)
    score_p <- logPosteriorBeta(beta=propose,y=YS,x=XS,zcond=zo,rho=rho[,i],u0S,V0invS,y_dist=select_dist)
    if (is.na(score_p-score_c)) {
      logA  <- -Inf
    } else {
      logA  <- min(0, score_p-score_c)
    }
    MH_acceptance_prob[i,cnames] <- logA
    if (logA > log(runif(1))) {
      betaS[,i] <- propose
      MH_acceptance[i,cnames] <- 1
    } else {
      betaS[,i] <- current
    }
    # print("update ZS")
    ## update zs: Gibbs from truncated normal ----
    ZS[,i] <- drawzDisc(y_dist=select_dist,y=YS,x=XS,beta=betaS[,i],zcond=zo,rho=rho[,i])
  }
  
  ##==== Return Iterates ====
  iterates <- list()
  iterates$rho <- rho[,-1]
  iterates$betaO <- betaO[,-1]
  iterates$betaS <- betaS[,-1]
  iterates$ZS <- ZS[,-1]
  if (outcome_dist=="Normal") {
    iterates$sigmaO <- sigmaO[,-1]
  } else if (outcome_dist=="Negative Binomial") {
    iterates$rO <- rO[,-1]
  } else {}
  iterates$MH_step <- MH_step
  iterates$MH_acceptance_prob <- MH_acceptance_prob
  
  return(iterates)
}
