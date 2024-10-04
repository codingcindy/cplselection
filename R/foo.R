# # library(usethis)
# # library(devtools)
# # load_all()
# 
# set.seed(1048)
# iterates <- cplselectionMCMC(outcome_formula, select_formula,
#                              outcome_dist, select_dist, data=df,
#                              loop, stepsize, stepadj, sliceadj)
# res <- cplselectionInfer(iterates, burnin)
# round(cbind(trueval, res), 2)