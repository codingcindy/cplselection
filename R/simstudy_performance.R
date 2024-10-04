#' Generate Performance Metrics Table for Simulation Studies
#'
#' @param restab Result summary table
#'
#' @importFrom utils head
#' @importFrom stats median 
#' @importFrom dplyr %>% group_by group_by_at mutate rename summarise ungroup  
#' 
#' @return A data.frame object containing performance metrics.
#' @export
#'
#' @examples NULL
simstudy_performance <- function(restab, by=c("x_dist","dependence","selectivity"), type="Bayes",trim=0) {

  ## Standardize Input Dataframe
  if (type=="Bayes") {
    ## add dependence & selectivity columns
    restab <- restab %>% 
      rename(estimate = median) %>%
      group_by(select_dist,outcome_dist,x_dist,parset,seed) %>% 
      mutate(dependence = trueval[parameter=="rho"], 
             selectivity = ifelse((x_dist=="Poisson"), 
                                  trueval[parameter=="betaS0"]+1,
                                  trueval[parameter=="betaS0"])) %>% 
      ungroup()
  } else if (type=="Standard") {
    ## add CI boundaries
    restab <- restab %>% 
      mutate(lb = estimate-1.96*se, ub = estimate+1.96*se)
  } else {
    stop("Type not supported.")
  }
  ## Generate Performance Metrics
  group <- c("select_dist","outcome_dist","parameter", by)
  perftab <- restab %>% 
    group_by_at(group) %>% 
    summarise(
      MAPE = mean(abs(estimate[trueval!=0]/trueval[trueval!=0]-1), na.rm=TRUE, trim=trim), 
      MAE = mean(abs(estimate-trueval), trim=trim), RMSE = sqrt(mean((estimate-trueval)^2, trim=trim)), 
      inCI = mean((trueval<=ub)&(trueval>=lb), trim=trim), widthCI = mean(ub-lb, trim=trim)
    )
  
  return(perftab)
}
