simstudy_ProbitNonnorm <- function(seeds, folderpath, rest.time=600) {
  library(parallel)
  ## create directory if not exists
  if (!file.exists(folderpath)) {
    dir.create(folderpath)
  }
  ## manipulations
  sd <- "Probit"
  od <- c("Probit", "Logit", "Poisson", "Negative Binomial")
  cor <- c(-.8,-.5,-.2,0,.2,.5,.8)
  selectivity <- c(-1,0,1) # low 30%, medium 50%, high 70%
  endog_dist_list <- list(list(dist="Normal",mean=0,sd=1),
                          list(dist="Uniform",min=-1,max=1),
                          list(dist="Poisson",lambda=1),
                          list(dist="Normal Mixture", weight=c(0.5,0.5), mean=c(-1,1), sd=c(1,2))
  )
  ## assignment
  nobs <- 500
  loop <- 20000
  burnin <- 10000
  stepsize <- 0.2
  stepadj <- 0.1
  sliceadj <- 100
  select_formula <- YS~X1+X2+X3
  outcome_formula <- YO~X1+X2
  outcome_par <- list(beta=c(1,1,1,0), sigma=1, r=1)
  parsetting <- 0
  time_start <- Sys.time()
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
            ## temporary use: add more outcome distributions
            message(paste0("Parameter setting ", parsetting, " starts."))
            if (outcome_dist != "Negative Binomial") next
            # ## temporary use: restart if stopped halfway
            # message(paste0("Parameter setting ", parsetting, " starts."))
            # if (parsetting < 135) next
            filename <- paste0("xdist",i,"vis",pct,"dep",theta,"_",select_dist,"-",outcome_dist)
            tryCatch(
              {
                parallel::mclapply(seeds, FUN=simcplEstimate, 
                                   filedir=folderpath, filename=filename,
                                   outcome_dist=outcome_dist, select_dist=select_dist, 
                                   outcome_par=outcome_par, select_par=select_par, 
                                   theta=theta, x_dist=x_dist, nobs=nobs, 
                                   outcome_formula=outcome_formula, select_formula=select_formula,
                                   loop=loop, burnin=burnin, stepsize=stepsize, stepadj=stepadj, sliceadj=sliceadj,
                                   mc.preschedule=TRUE, mc.cores=4)
              }, 
              error = function(cond) {
                message(conditionMessage(cond))
                NA
              },
              finally = {
                message(paste0("Parameter setting ",parsetting," done."))
                # ## avoid overheating, silence if run on server
                # time_run <- Sys.time()-time_start
                # if (time_run>3000) {
                #   message("Running for ",time_run/60, " minutes. Take a rest for ", rest.time/60," minutes.")
                #   Sys.sleep(rest.time)
                #   time_start <- Sys.time()
                # }
              }
            )
          }
        }
      }
    }
  }
  return(NULL)
}