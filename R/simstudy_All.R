#' Run all simulation studies with cplselection functions
#'
#' @param seeds random seeds 
#' @param folderpath directory to save the .RData results
#' @param rest.time if to avoid overheating on personal device, time of resting 
#' @param mc.cores parallel cores
#'
#' @return NULL
#' @export
#'
#' @examples NULL
simstudy_All <- function(seeds, folderpath, rest.time=0, mc.cores=8) {
  library(parallel)
  ## create directory if not exists
  if (!file.exists(folderpath)) {
    dir.create(folderpath)
  }
  ## manipulations
  sd <- c("Probit", "Logit")
  od <- c("Normal", "Probit", "Logit", "Poisson", "Negative Binomial")
  # od <- c("Exoponential")
  cor <- c(-.8,-.5,-.2,0,.2,.5,.8)
  vis <- c(0.2, 0.5, 0.8) # low, medium, high
  nobs <- c(100, 500) # low, high
  endog_dist_list <- list(list(dist="Normal",mean=0,sd=1),
                          list(dist="Uniform",min=-1,max=1),
                          list(dist="Poisson",lambda=1),
                          list(dist="Normal Mixture", weight=c(0.5,0.5), mean=c(-1,1), sd=c(1,2))
  )
  ## assignment
  loop <- 20000
  burnin <- 10000
  stepsize <- 0.2
  stepadj <- 0.1
  sliceadj <- 100
  select_formula <- YS~X1+X2+X3
  outcome_formula <- YO~X1+X2
  outcome_par <- list(beta=c(1,1,1,0), sigma=1)
  select_par <- list(beta=c(NA,1,0,1))
  parsetting <- 0
  for (select_dist in sd) {
    for (outcome_dist in od) {
      for (n in nobs) {
        for (v in vis) {
          for (theta in cor) {
            for (i in (1:length(endog_dist_list))) {
              x_dist <- list(endog_dist_list[[i]],
                             list(dist="Normal",mean=0,sd=1),
                             list(dist="Normal",mean=0,sd=1))
              parsetting <- parsetting+1
              filename <- paste0(select_dist,outcome_dist,"_N",n,"Vis",100*v,"Dep",theta*10,"_X1",x_dist[[1]]$dist)
              tryCatch(
                {
                  parallel::mclapply(seeds, FUN=simcplEstimate,
                                     filedir=folderpath, filename=filename,
                                     outcome_dist=outcome_dist, select_dist=select_dist,
                                     outcome_par=outcome_par, select_par=select_par,
                                     theta=theta, x_dist=x_dist, nobs=n, vis=v,
                                     outcome_formula=outcome_formula, select_formula=select_formula,
                                     loop=loop, burnin=burnin, stepsize=stepsize, stepadj=stepadj, sliceadj=sliceadj,
                                     mc.preschedule=TRUE, mc.cores=mc.cores)
                }, 
                error = function(cond) {
                  message(conditionMessage(cond))
                  NA
                },
                finally = {
                  # ## if overheating is a concern
                  # message(paste0("Parameter setting ",parsetting," done. Rest for ",rest.time/60," minutes."))
                  # Sys.sleep(rest.time)
                  message(paste0("Parameter setting ",parsetting," done."))
                }
              )
            }
          }
        }
      }
    }
  }
  return(NULL)
}