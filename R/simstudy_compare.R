# simstudy_compare <- function(sumtab_heckit, sumtab_cpl, parname, by, trim=0) {
#   ## summarize performance
#   # perftab_ols <- cbind(method="OLS", simstudy_performance(sumtab_ols, type="Standard", by=by, trim=trim)
#   perftab_heckit <- cbind(method="Heckit", simstudy_performance(sumtab_heckit,type="Standard", by=by, trim=trim))
#   perftab_cpl <- cbind(method="Copula", simstudy_performance(sumtab_cpl, by=by, trim=trim))
#   ## select parameter, order columns, sort
#   colorder <- c("select_dist","outcome_dist","method",by,"MAE","RMSE","inCI","widthCI")
#   perftab <- rbind(
#     # perftab_ols[perftab_ols$parameter==parname, colorder],
#     perftab_heckit[perftab_heckit$parameter==parname, colorder], 
#     perftab_cpl[perftab_cpl$parameter==parname, colorder])
#   perftab <- arrange(perftab, desc(select_dist), outcome_dist, desc(method))
#   
#   return(perftab)
# }