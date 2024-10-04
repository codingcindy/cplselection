#' Summarize the Bayesian estimation results of input files
#'
#' @importFrom utils head
#' @return A summary table of Bayesian estimates and true parameters
#' @export 
#'
#' @examples NULL
simstudy_summary <- function(folderpath) {
  ## get file names
  file_list <- list.files(folderpath)
  ## load files
  sumtab <- data.frame()
  for (f in file_list) {
    load(paste0(folderpath,"/",f))
    if (reslist$outcome_dist=="Normal"|reslist$outcome_dist=="Negative Binomial") {
      trueval <- c(reslist$theta, head(reslist$outcome_par$beta, -1), 
                   reslist$select_par$beta,reslist$outcome_par[[2]])
    } else {
      trueval <- c(reslist$theta, head(reslist$outcome_par$beta, -1), 
                   reslist$select_par$beta)
    }
    tmptab <- data.frame(select_dist=reslist$select_dist,
                         outcome_dist=reslist$outcome_dist,
                         x_dist=reslist$x_dist[[1]]$dist,
                         parameter=rownames(reslist$result), 
                         trueval=trueval, reslist$result, 
                         fname=f, row.names=NULL)
    sumtab <- rbind(sumtab, tmptab)
  }
  tmpstr <- do.call(rbind, strsplit(gsub("parset|seed|\\.RData", "", sumtab$fname), "_"))
  sumtab$parset <- tmpstr[,1]
  sumtab$seed <- tmpstr[,2]
  ## reorder columns
  cnames <- c("select_dist","outcome_dist","x_dist","parset","seed",
              "parameter","trueval","median","sd","lb","ub")
  sumtab <- sumtab[,cnames]
  
  return(sumtab)
}