# ## Testing Self-defined R Functions
# rm(list=ls())
# library(usethis)
# library(devtools)
# library(here)
# load_all()
# 
# ####------ Inputs ------
# ## marginal specifications
# outcome_dist <- "Negative Binomial"
# select_dist <- "Probit"
# outcome_par <- list(beta=c(1,1,1,0), sigma=1, r=1)
# select_par <- list(beta=c(NA,1,0,1))
# ## key manipulations
# theta <- -0.8
# nobs <- 100
# vis <- 0.25
# x_dist <- list(list(dist="Uniform",min=-1,max=1),
#                list(dist="Normal",mean=0,sd=1),
#                list(dist="Normal",mean=0,sd=1))
# # ## X1 alternatives:
# #   list(dist="Uniform",min=-1,max=1)
# #   list(dist="Normal",mean=0,sd=1)
# #   list(dist="Poisson",lambda=1)
# #   list(dist="Normal Mixture", weight=c(0.5,0.5), mean=c(-1,1), sd=c(1,2))
# #   list(dist="Negative Binomial", mu=exp(0), size=1)
# 
# ####------ Test: simdata ------
# seed <- 2022
# set.seed(seed)
# sim <- simdata(outcome_dist=outcome_dist, outcome_par=outcome_par,
#                select_dist=select_dist, select_par=select_par,
#                theta=theta, nobs=nobs, vis=vis, x_dist=x_dist,
#                selection=TRUE)
# df <- sim$data
# select_par$beta[1] <- sim$intercept
# ## checks
# (mean(sim$data$YS)==vis)
# summary(df)
# 
# #### ------ Test: Estimation ------
# #### cplselectionMCMC, cplselectionInfer, simHeckitEstimate
# ## Estimate
# set.seed(seed)
# iterates <- cplselectionMCMC(outcome_formula, select_formula,
#                              outcome_dist, select_dist, data=df,
#                              loop, stepsize, stepadj, sliceadj)
# ## Infer
# if (outcome_dist=="Normal") {
#   # morepar <- outcome_par$sigma
#   trueval <- c(rho=theta, betaO=outcome_par$beta[1:3],
#                betaS=select_par$beta, sigma=outcome_par$sigma)
# } else if (outcome_dist=="Negative Binomial") {
#   # morepar <- outcome_par$r
#   trueval <- c(rho=theta, betaO=outcome_par$beta[1:3],
#                betaS=select_par$beta, r=outcome_par$r)
# } else {
#   trueval <- c(rho=theta, betaO=outcome_par$beta[1:3],
#                betaS=select_par$beta)
# }
# res <- cplselectionInfer(iterates, burnin, trueval)
# res_cpl <- round(res, 2)
# 
# ## Compare
# # heck style
# heck <- simHeckitEstimate(
#   seed=seed, outcome_dist=outcome_dist, select_dist=select_dist,
#   outcome_par=outcome_par, select_par=select_par, theta=theta,
#   x_dist=x_dist, nobs=nobs, vis=vis,
#   outcome_formula=outcome_formula, select_formula=select_formula)
# res_heck <- round(heck[c(5:7,1:4,9),c("trueval","estimate","se","pval")], 3)
# # naive regression
# if (outcome_dist=="Probit") {
#   naive <- glm(outcome_formula, data=df, family=binomial(link="probit"))
# } else if (outcome_dist=="Logit") {
#   naive <- glm(outcome_formula, data=df, family=binomial(link="probit"))
# } else if (outcome_dist=="C-log-log") {
#   naive <- glm(outcome_formula, data=df, family=poisson(link="log"))
# } else if (outcome_dist=="Normal") {
#   naive <- lm(outcome_formula, data=df)
# } else if (outcome_dist=="Poisson") {
#   naive <- MASS::glm.nb(outcome_formula, data=df)
# } else if (outcome_dist=="Negative Binomial") {
#   naive <- MASS::glm.nb(outcome_formula, data=df)
# } else {}
# res_naive <- c(naive$coefficients,NA,NA,NA,NA)
# names(naive) <- c("betoO0","betaO1","betaO2","betoS0","betaS1","betaS2","betaS3")
# 
# # tablist
# round(cbind(res_cpl, naive=c(NA,res_naive,NA),
#             heck2s=heck$estimate[c(10,5:7,1:4,9)]), 2)
# 
# ####------ Test: Bulk Estimation ------
# #### simcplEstimate
# ## (i) save reslist under folder tmpdir and
# ## (ii) return a result summary to tmp
# tmpdir <- paste0(here(), "/data")
# tmp <- simcplEstimate(
#   seed=1, filedir=tmpdir, filename="delete",
#   outcome_dist, select_dist, outcome_par, select_par, theta, x_dist, nobs, vis,
#   outcome_formula, select_formula, loop, burnin, stepsize, stepadj, sliceadj
# )
# round(tmp, 2)
# 
# #### simstudy_ProbitNorm, simstudy_ProbitNonnorm, simstudy_LogitNorm, simstudy_LogitNonnorm, simstudy_All
# ## (i) save RData files under folderpath
# ## (ii) return NULL
# seeds <- 1:2
# rest.time <- 0
# mc.cores <- 8
# folderpath <- paste0(here(), "/data/simstudy_ProbitNorm")
# # tests
# simstudy_ProbitNorm(
#   seeds=seeds, folderpath=folderpath, rest.time=rest.time, mc.cores=mc.cores)
# 
# simstudy_All(
#   seeds=seeds, folderpath=folderpath, rest.time=rest.time, mc.cores=mc.cores)
# 
