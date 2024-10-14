# ## This file performs main analysis of simulation studies
# rm(list=ls())
# # library(cplselection)
# library(usethis)
# library(devtools)
# load_all()
# 
# ## Study 1:
# ##    step 1. simulate data set
# ##    ** DGP parameters
# ##    ** b0: nonzero intercept
# ##    ** x1: endogenous var
# ##    ** x2: exogenous var
# ##    ** x3: IV-type var
# 
# ## Works well:
# ##    Normal - Probit/Logit/Cloglog
# ##    Poisson - Probit/Logit/Cloglog
# ## Not working:
# ##    NB - Probit
# outcome_dist <- "Negative Binomial"
# select_dist <- "Probit"
# outcome_par <- list(beta=c(1,1,1,0), sigma=1, r=1)
# select_par <- list(beta=c(-1,1,0,1)) #
# theta <- -0.8
# nobs <- 500
# seed <- 2022
# x_dist <- list(list(dist="Poisson",lambda=1),
#                list(dist="Normal",mean=0,sd=1),
#                list(dist="Normal",mean=0,sd=1))
# # x_dist <- list(list(dist="Normal",mean=0,sd=1),
# #                list(dist="Uniform",min=-1,max=1),
# #                list(dist="Poisson",lambda=1),
# #                list(dist="Normal Mixture",
# #                     weight=c(0.5,0.5), mean=c(-1,1), sd=c(1,2))
# #                # list(dist="Negative Binomial", size=5, prob=0.8)
# #                )
# if (outcome_dist=="Normal") {
#   morepar <- outcome_par$sigma
# } else if (outcome_dist=="Negative Binomial") {
#   morepar <- outcome_par$r
# } else {
#   morepar <- NULL
# }
# trueval <- c(theta, outcome_par$beta[1:3], select_par$beta, morepar)
# 
# ##    ** Generate Dataset
# set.seed(seed)
# df <- simdata(outcome_dist=outcome_dist, outcome_par=outcome_par,
#               select_dist=select_dist, select_par=select_par,
#               theta=theta, x_dist=x_dist, nobs=nobs, selection=TRUE)
# summary(df$YS)
# ##    step 2. copula MCMC iterations
# ##    ** specify margins
# outcome_formula <- YO~X1+X2
# select_formula <- YS~X1+X2+X3
# ##    ** MCMC settings
# loop <- 20000
# burnin <- 5000
# stepsize <- 0.2
# stepadj <- 0.1
# sliceadj <- 100
# 
# set.seed(seed)
# iterates <- cplselectionMCMC(outcome_formula, select_formula,
#                              outcome_dist, select_dist, data=df,
#                              loop, stepsize, stepadj, sliceadj)
# res <- cplselectionInfer(iterates, burnin)
# round(cbind(trueval, res), 2)
# ## why not updating in NB case ? logA always -Inf
# ## logposteriorbeta not revised in nonnormal cases!
# 
# 
# # compare to binary regression
# binom <- glm(select_formula, data=df,
#              family=binomial(link="probit"))
# nb <- MASS::glm.nb(outcome_formula, data = df)
# # naive <- lm(outcome_formula, data=df)
# # naive <- glm(outcome_formula, data=df,
# #              family=poisson(link="log"))
# naive <- cbind(nbreg=c(nb$coefficients,NA,NA,NA,NA),
#                probit=c(NA,NA,NA,binom$coefficients))
# rownames(naive) <- c("betoO0","betaO1","betaO2",
#                      "betoS0","betaS1","betaS2","betaS3")
# 
# heck <- simHeckitEstimate(
#   seed=seed, outcome_dist, select_dist,
#   outcome_par, select_par, theta, x_dist, nobs,
#   outcome_formula, select_formula)
# round(heck[c(5:7,1:4,9),c("trueval","estimate","se","pval")], 3)
# 
# round(cbind(trueval, res, rbind(NA,naive,NA),
#             heck2s=heck$estimate[c(10,5:7,1:4,9)]), 2)
# 
# #
# # # Diagnostics
# # par(mfrow=c(3,3))
# # plotIters(iterates,"betaS")
# # plotIters(iterates,"betaO")
# # plotIters(iterates,"rho")
# # plotIters(iterates,"sigmaO")
# # plotIters(iterates,"rO")
# #
# # par(mfrow=c(2,3))
# # plot(iterates$MH_step[,"rho"], type="l",main="rho")
# # plot(iterates$MH_step[,"betaS0"], type="l",main="betaS0")
# # plot(iterates$MH_step[,"betaS1"], type="l",main="betaS1")
# # plot(iterates$MH_step[,"betaS2"], type="l",main="betaS2")
# # plot(iterates$MH_step[,"betaS3"], type="l",main="betaS3")
# # plot(iterates$MH_step[,"betaO0"], type="l",main="betaO0")
# # plot(iterates$MH_step[,"betaO1"], type="l",main="betaO1")
# # plot(iterates$MH_step[,"betaO2"], type="l",main="betaO2")
# # plot(iterates$MH_step[,"sigmaO"], type="l",main="sigmaO")
# # plot(iterates$MH_step[,"rO"], type="l",main="rO")
# #
# # par(mfrow=c(2,3))
# # plot(iterates$MH_acceptance_prob[,"rho"], type="l",main="rho")
# # plot(iterates$MH_acceptance_prob[,"betaS0"], type="l",main="betaS0")
# # plot(iterates$MH_acceptance_prob[,"betaS1"], type="l",main="betaS1")
# # plot(iterates$MH_acceptance_prob[,"betaS2"], type="l",main="betaS2")
# # plot(iterates$MH_acceptance_prob[,"betaS3"], type="l",main="betaS3")
# # plot(iterates$MH_acceptance_prob[,"betaO0"], type="l",main="betaO0")
# # plot(iterates$MH_acceptance_prob[,"betaO1"], type="l",main="betaO1")
# # plot(iterates$MH_acceptance_prob[,"betaO2"], type="l",main="betaO2")
# # plot(iterates$MH_acceptance_prob[,"sigmaO"], type="l",main="sigmaO")
# # plot(iterates$MH_acceptance_prob[,"rO"], type="l",main="rO")
# #
# # ## Mass Estimation
# # rm(list=ls())
# # library(usethis)
# # library(devtools)
# # load_all()
# #
# # seeds <- 1:100
# # folderpath <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudyProbitNorm"
# # simstudy_ProbitNorm(seeds, folderpath, rest.time=600)
# #
# # folderpath <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudyProbitNonnorm"
# # simstudy_ProbitNonnorm(seeds, folderpath, rest.time=360)
# #
# # folderpath <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudyLogitNnorm"
# # simstudy_LogitNorm(seeds, folderpath, rest.time=180)
# #
# # folderpath <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudyLogitNonnorm"
# # simstudy_LogitNonnorm(seeds, folderpath, rest.time=180)
# #
# #
# # ## Performance Summary
# # rm(list=ls())
# # library(devtools)
# # library(roxygen2)
# # library(usethis)
# # load_all()
# # #  load and summarize data (takes ~10 minutes)
# # folderpath <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudyProbitNorm"
# # sumtab_probitnorm <- simstudy_summary(folderpath)
# #
# # folderpath <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudyProbitNonnorm"
# # sumtab_probitnonnorm <- simstudy_summary(folderpath)
# #
# # folderpath <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudyLogitNnorm"
# # sumtab_logitnorm <- simstudy_summary(folderpath)
# #
# # folderpath <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudyLogitNonnorm"
# # sumtab_logitnonnorm <- simstudy_summary(folderpath)
# #
# # #  save to
# # prefix <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudy"
# # save(sumtab_probitnorm, sumtab_probitnonnorm, sumtab_logitnorm, sumtab_logitnonnorm,
# #      file = paste0(prefix, "/sumtab_bayes.RData"))
# #
# # # overall table: no trim
# # prefix <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudy"
# # load(paste0(prefix, "/sumtab_bayes.RData"))
# # write.table(simstudy_performance(sumtab_probitnorm, by=c()),
# #             sep=",", row.names=FALSE, col.names=TRUE, append=FALSE,
# #             file=paste0(prefix,"/performance_bayes.csv"))
# # write.table(simstudy_performance(sumtab_probitnonnorm, by=c()),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes.csv"))
# # write.table(simstudy_performance(sumtab_logitnorm, by=c()),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes.csv"))
# # write.table(simstudy_performance(sumtab_logitnonnorm, by=c()),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes.csv"))
# # # detail tables: no trim
# # write.table(simstudy_performance(sumtab_probitnorm),
# #             sep=",", row.names=FALSE, col.names=TRUE, append=FALSE,
# #             file=paste0(prefix,"/performance_bayes_detail.csv"))
# # write.table(simstudy_performance(sumtab_probitnonnorm),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes_detail.csv"))
# # write.table(simstudy_performance(sumtab_logitnorm),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes_detail.csv"))
# # write.table(simstudy_performance(sumtab_logitnonnorm),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes_detail.csv"))
# #
# # ## bayes performance: trimmed
# # trim <- 0.01
# # prefix <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudy"
# # load(paste0(prefix, "/sumtab_bayes.RData"))
# # write.table(simstudy_performance(sumtab_probitnorm, by=c(), trim=trim),
# #             sep=",", row.names=FALSE, col.names=TRUE, append=FALSE,
# #             file=paste0(prefix,"/performance_bayes_trim.csv"))
# # write.table(simstudy_performance(sumtab_probitnonnorm, by=c(), trim=trim),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes_trim.csv"))
# # write.table(simstudy_performance(sumtab_logitnorm, by=c(), trim=trim),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes_trim.csv"))
# # write.table(simstudy_performance(sumtab_logitnonnorm, by=c(), trim=trim),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes_trim.csv"))
# # # detail tables: trimmed
# # write.table(simstudy_performance(sumtab_probitnorm, trim=trim),
# #             sep=",", row.names=FALSE, col.names=TRUE, append=FALSE,
# #             file=paste0(prefix,"/performance_bayes_detail_trim.csv"))
# # write.table(simstudy_performance(sumtab_probitnonnorm, trim=trim),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes_detail_trim.csv"))
# # write.table(simstudy_performance(sumtab_logitnorm, trim=trim),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes_detail_trim.csv"))
# # write.table(simstudy_performance(sumtab_logitnonnorm, trim=trim),
# #             sep=",", row.names=FALSE, col.names=FALSE, append=TRUE,
# #             file=paste0(prefix,"/performance_bayes_detail_trim.csv"))
# #
# # ## Simulation Study with Heckit
# # rm(list=ls())
# # library(devtools)
# # library(roxygen2)
# # library(usethis)
# # load_all()
# # # ~15 mins, avoid repetition
# # folderpath <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudyHeckit"
# # prefix <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudy"
# # # seeds <- 1:100
# # # system.time(sumtab_heckit <- simstudy_Heckit(seeds, folderpath, rest.time=0))
# # # save(sumtab_heckit, file=paste0(prefix, "/sumtab_heckit.RData"))
# # #  load
# # load(file=paste0(prefix, "/sumtab_heckit.RData"))
# # #  performance tables: no trim
# # write.table(simstudy_performance(sumtab_heckit,type="Standard",by=c()),
# #             sep=",", row.names=FALSE, col.names=TRUE,
# #             file=paste0(prefix,"/performance_heckit.csv"))
# # write.table(simstudy_performance(sumtab_heckit,type="Standard"),
# #             sep=",", row.names=FALSE, col.names=TRUE,
# #             file=paste0(prefix,"/performance_heckit_detail.csv"))
# # #  performance tables: trim
# # trim <- 0.01
# # write.table(simstudy_performance(sumtab_heckit,type="Standard",by=c(),trim=trim),
# #             sep=",", row.names=FALSE, col.names=TRUE,
# #             file=paste0(prefix,"/performance_heckit_trim.csv"))
# # write.table(simstudy_performance(sumtab_heckit,type="Standard",trim=trim),
# #             sep=",", row.names=FALSE, col.names=TRUE,
# #             file=paste0(prefix,"/performance_heckit_detail_trim.csv"))
# #
# # ## Simulation Study with OLS : Benchmark
# # library(devtools)
# # library(roxygen2)
# # library(usethis)
# # load_all()
# # prefix <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudy"
# # # takes ~100 seconds
# # seeds <- 1:100
# # system.time(sumtab_ols <- simstudy_OLS(seeds, prefix, rest.time=0))
# # save(sumtab_ols, file=paste0(prefix, "/sumtab_ols.RData"))
# #
# # ## ---- Check performance by manipulation ---------------------
# # rm(list=ls())
# # library(dplyr)
# # library(openxlsx)
# # library(devtools)
# # library(roxygen2)
# # library(usethis)
# # load_all()
# # ## ---- load data
# # trim <- 0.0  # trimming is a bad practice: nobody implements as such. do not use!
# # prefix <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudy"
# # load(file=paste0(prefix, "/sumtab_heckit.RData"))
# # load(file=paste0(prefix, "/sumtab_bayes.RData"))
# # load(file=paste0(prefix, "/sumtab_ols.RData"))
# # sumtab_cpl <- rbind(sumtab_probitnorm,sumtab_probitnonnorm,
# #                     sumtab_logitnorm,sumtab_logitnonnorm)
# # ## ---- Overall:
# # perftab <- simstudy_compare(
# #   sumtab_heckit, sumtab_cpl, by=c(), parname="betaO1", trim=trim) %>%
# #   arrange(desc(select_dist),outcome_dist,desc(method))
# # ## ---- By Dependence:
# # perftab_rho <- simstudy_compare(
# #   sumtab_heckit, sumtab_cpl, by="dependence", parname="betaO1", trim=trim) %>%
# #   arrange(method,dependence,select_dist,outcome_dist)
# # ## ---- By Selectivity:
# # perftab_s <- simstudy_compare(
# #   sumtab_heckit, sumtab_cpl, by="selectivity", parname="betaO1", trim=trim) %>%
# #   arrange(method,selectivity,select_dist,outcome_dist)
# # ## ---- By X1 distribution:
# # perftab_xdist <- simstudy_compare(
# #   sumtab_heckit, sumtab_cpl, by="x_dist", parname="betaO1", trim=trim) %>%
# #   arrange(method,x_dist,select_dist,outcome_dist)
# # ## ---- By Dependence X Selectivity X X1_dist
# # perftab_all <- simstudy_compare(
# #   sumtab_heckit, sumtab_cpl, by=c("x_dist","dependence","selectivity"),
# #   parname="betaO1", trim=trim)
# # ## --- Write to XLSX
# # wb <- createWorkbook()
# # addWorksheet(wb, "overall")
# # addWorksheet(wb, "by_dependence")
# # addWorksheet(wb, "by_selectivity")
# # addWorksheet(wb, "by_x1_dist")
# # addWorksheet(wb, "by_all_dimension")
# # writeData(wb, sheet="overall", x=perftab)
# # writeData(wb, sheet="by_dependence", x=perftab_rho)
# # writeData(wb, sheet="by_selectivity", x=perftab_s)
# # writeData(wb, sheet="by_x1_dist", x=perftab_xdist)
# # writeData(wb, sheet="by_all_dimension", x=perftab_all)
# # saveWorkbook(
# #   wb, file=paste0(prefix, "/perftab_manipulation.xlsx"), overwrite=TRUE)
# # # saveWorkbook(
# # #   wb, file=paste0(prefix, "/perftab_manipulation_trim001.xlsx"), overwrite=TRUE)
# # # saveWorkbook(
# # #   wb, file=paste0(prefix, "/perftab_manipulation_trim002.xlsx"), overwrite=TRUE)
# # # saveWorkbook(
# # #   wb, file=paste0(prefix, "/perftab_manipulation_trim010.xlsx"), overwrite=TRUE)
# #
# # ## ---- Relative Performance: t-test
# # df_test <- perftab_all %>%
# #   tidyr::pivot_wider(names_from=method, values_from=c(MAE,RMSE,inCI,widthCI))
# #
# # res_test <- df_test %>%
# #   group_by(select_dist, outcome_dist) %>%
# #   summarise(
# #     MAE_est = t.test(MAE_Heckit, MAE_Copula, paired=TRUE)$estimate,
# #     MAE_t = t.test(MAE_Heckit, MAE_Copula, paired=TRUE)$statistic,
# #     MAE_p = t.test(MAE_Heckit, MAE_Copula, paired=TRUE)$p.value,
# #     RMSE_est = t.test(RMSE_Heckit, RMSE_Copula, paired=TRUE)$estimate,
# #     RMSE_t = t.test(RMSE_Heckit, RMSE_Copula, paired=TRUE)$statistic,
# #     RMSE_p = t.test(RMSE_Heckit, RMSE_Copula, paired=TRUE)$p.value,
# #     inCI_est = t.test(inCI_Heckit, inCI_Copula, paired=TRUE)$estimate,
# #     inCI_t = t.test(inCI_Heckit, inCI_Copula, paired=TRUE)$statistic,
# #     inCI_p = t.test(inCI_Heckit, inCI_Copula, paired=TRUE)$p.value,
# #     widthCI_est = t.test(widthCI_Heckit, widthCI_Copula, paired=TRUE)$estimate,
# #     widthCI_t = t.test(widthCI_Heckit, widthCI_Copula, paired=TRUE)$statistic,
# #     widthCI_p = t.test(widthCI_Heckit, widthCI_Copula, paired=TRUE)$p.value
# #   ) %>%
# #   mutate(across(where(is.numeric), round, 4))
# #
# # ## ---- Regression / ANOVA Analysis of Bias Factors ----------
# # rm(list=ls())
# # library(usethis)
# # library(devtools)
# # library(roxygen2)
# # load_all()
# # library(dplyr)
# # library(broom)
# # # ## ---- prepare regression dataset
# # # prefix <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudy"
# # # load(file=paste0(prefix, "/sumtab_OLS.RData"))
# # # load(file=paste0(prefix, "/sumtab_heckit.RData"))
# # # load(file=paste0(prefix, "/sumtab_bayes.RData"))
# # # dfmerge <- sumtab_heckit %>%
# # #   dplyr::select(select_dist,outcome_dist,x_dist,dependence,selectivity,seed,pctobs) %>%
# # #   unique()
# # # dfols <- sumtab_ols %>%
# # #   dplyr::mutate(method="OLS", lb=estimate-1.96*se, ub=estimate+1.96*se)
# # # dfheck <- sumtab_heckit %>%
# # #   dplyr::mutate(method="Heckit", lb=estimate-1.96*se, ub=estimate+1.96*se)
# # # dfreg <- rbind(sumtab_probitnorm,sumtab_probitnonnorm,
# # #                sumtab_logitnorm,sumtab_logitnonnorm) %>%
# # #   dplyr::rename(estimate=median, se=sd) %>%
# # #   dplyr::group_by(select_dist,outcome_dist,x_dist,seed,parset) %>%
# # #   dplyr::mutate(method="Copula", seed=as.numeric(seed),
# # #                 dependence=trueval[parameter=="rho"],
# # #                 selectivity=ifelse(x_dist=="Poisson",
# # #                                    trueval[parameter=="betaS0"]+1,
# # #                                    trueval[parameter=="betaS0"]
# # #                 )) %>%
# # #   dplyr::ungroup() %>% dplyr::select(-parset) %>%
# # #   dplyr::left_join(dfmerge, by=c("select_dist","outcome_dist","x_dist","dependence","selectivity","seed")) %>%
# # #   dplyr::bind_rows(dfols) %>%
# # #   dplyr::bind_rows(dfheck) %>%
# # #   dplyr::mutate( # performance metrics
# # #     AE = abs(estimate - trueval),
# # #     SE = AE^2,
# # #     inCI = (trueval>=lb)*(trueval<=ub),
# # #     lenCI = (ub - lb)
# # #   )
# # # rm(list=ls()[grep(pattern="sumtab*", x=ls())])
# # # rm(dfmerge, dfheck, dfols)
# # # save(dfreg, file=paste0(prefix,"/bias_factor_df.RData"))
# # ## ---- run regression for impact factors on endog var coef
# # library(modelsummary)
# # prefix <- "/Users/xindieh/Documents/3_Academics/Projects/8_SelectionBias/simulation/simstudy"
# # load(file=paste0(prefix,"/bias_factor_df.RData"))
# # filter <- 0.003 ## remove outliers for each method
# # df <- dfreg %>%
# #   dplyr::filter(parameter=="betaO1") %>%
# #   dplyr::group_by(method) %>%
# #   dplyr::filter(AE<=quantile(AE, 1-filter)) %>%
# #   dplyr::ungroup()
# # mods <- df %>%
# #   dplyr::group_by(select_dist, outcome_dist, method) %>%
# #   dplyr::do(
# #     AE = lm(AE ~ (abs(dependence) + selectivity + factor(x_dist)), data = .),
# #     SE = lm(SE ~ (abs(dependence) + selectivity + factor(x_dist)), data = .),
# #     inCI = glm(inCI ~ abs(dependence) + selectivity + factor(x_dist), family = binomial(link="logit"), data = .),
# #     lenCI = lm(lenCI ~ abs(dependence) + selectivity + factor(x_dist), data = .)
# #   )
# # coef_map <- c("abs(dependence)"="Dependence",
# #               "selectivity"="Visibility",
# #               "factor(x_dist)Uniform"="~ uniform",
# #               "factor(x_dist)Poisson"="~ poisson",
# #               "factor(x_dist)Normal Mixture"="~ normal mixture")
# # gof_map <- c("nobs", "adj.r.squared")
# # rows <- tribble(~term,~Copula,~Heckit,~OLS,~Copula,~Heckit,~OLS,~Copula,~Heckit,~OLS,~Copula,~Heckit,~OLS,
# #                 "X1 ~ normal","-","-","-","-","-","-","-","-","-","-","-","-")
# # attr(rows,"position") <- 7
# # tables <- list()
# # nmods <- 3
# # for (m in (1:nrow(mods))) {
# #   if (m%%nmods==1) {
# #     rangemods <- c(m:(m+nmods-1))
# #     outlist <- list(
# #       "Absolute Error" = unlist(mods[rangemods,"AE"], recursive=FALSE) %>% setNames(mods$method[rangemods]),
# #       "Squared Error" = unlist(mods[rangemods,"SE"], recursive=FALSE) %>% setNames(mods$method[rangemods]),
# #       "True Value in 95%CI" = unlist(mods[rangemods,"inCI"], recursive=FALSE) %>% setNames(mods$method[rangemods]),
# #       "Length of 95%CI" = unlist(mods[rangemods,"lenCI"], recursive=FALSE) %>% setNames(mods$method[rangemods])
# #     )
# #     main <- paste0("Selection - Outcome : ", mods$select_dist[m], " - ", mods$select_dist[m])
# #     tables[[1+m%/%nmods]] <- modelsummary(
# #       outlist, stars=TRUE, title=main, add_rows=rows,
# #       align="ldddddddddddd", shape="cbind", coef_map=coef_map, gof_map=gof_map)
# #   } else {}
# # }
# # ## ---- anova plots
# # library(ggpubr)
# # library(viridis)
# #
# #
