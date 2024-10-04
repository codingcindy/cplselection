#' Performs simulation study with heckit method
#'
#' @param seeds 
#' @param folderpath 
#' @param rest.time 
#'
#' @return A data.frame object containing results.
#' @export
#'
#' @examples NULL
simstudy_OLS <- function(seeds, folderpath, rest.time) {
  library(parallel)
  ## create directory if not exists
  if (!file.exists(folderpath)) {
    dir.create(folderpath)
  }
  ## manipulations
  sd <- c("Probit", "Logit")
  od <- c("Normal", "Probit", "Logit", "Poisson")
  cor <- c(-.8,-.5,-.2,0,.2,.5,.8)
  selectivity <- c(-1,0,1) # low 30%, medium 50%, high 70%
  endog_dist_list <- list(list(dist="Normal",mean=0,sd=1),
                          list(dist="Uniform",min=-1,max=1),
                          list(dist="Poisson",lambda=1),
                          list(dist="Normal Mixture", weight=c(0.5,0.5), mean=c(-1,1), sd=c(1,2))
  )
  ## assignment
  nobs <- 500
  select_formula <- YS~X1+X2+X3
  outcome_formula <- YO~X1+X2
  outcome_par <- list(beta=c(1,1,1,0), sigma=1)
  t0 <- Sys.time()
  parsetting <- 0
  resdf <- data.frame()
  for (i in (1:length(endog_dist_list))) {
    x_dist <- list(endog_dist_list[[i]],
                   list(dist="Normal",mean=0,sd=1),
                   list(dist="Normal",mean=0,sd=1))
    for (pct in selectivity) {
      if (endog_dist_list[[i]]$dist=="Poisson") {
        select_par <- list(beta=c(pct-1,1,0,1))
      } else {
        select_par <- list(beta=c(pct,1,0,1))
      }
      for (select_dist in sd) {
        for (outcome_dist in od) {
          for (theta in cor) {
            parsetting <- parsetting+1
            tryCatch(
              {
                reslist <- parallel::mclapply(
                  seeds, FUN=simOLSEstimate, 
                  outcome_dist=outcome_dist, select_dist=select_dist, 
                  outcome_par=outcome_par, select_par=select_par, 
                  theta=theta, x_dist=x_dist, nobs=nobs, 
                  outcome_formula=outcome_formula, select_formula=select_formula,
                  mc.preschedule=TRUE, mc.cores=8)
              }, 
              error = function(cond) {
                message(conditionMessage(cond))
                NA
              },
              finally = {
                # save(reslist, file=paste0(folderpath, "/parset", parsetting, ".RData"))
                tmpdf <- do.call(rbind, reslist) 
                resdf <- rbind(resdf, tmpdf)
                if (as.numeric(Sys.time()-t0)%/%1200>0) {
                  t0 <- Sys.time()
                  message(paste0("Running over 20 mins. Rest for ",rest.time/60," minutes."))
                  Sys.sleep(rest.time)
                }
              }
            )
          }
        }
      }
    }
  }
  ## automatic save
  save(resdf, file=paste0(folderpath, "/simstudy_OLS.RData"))
  return(resdf)
}