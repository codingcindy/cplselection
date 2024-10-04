plotIters <- function(iterates, parnames=NULL) {
  g <- list()
  if (is.null(parnames)) {
    parnames <- parnames[parnames!="ZS"]
  } else {
    for (par in parnames) {
      if (is.null(dim(iterates[[par]]))) {
        g[[par]] <- plot(iterates[[par]], type="l", xlab="iteration", main=par)
      } else {
        for (k in (1:(dim(iterates[[par]])[1]))) {
          yname <- paste0(par,k-1)
          g[[yname]] <- plot(iterates[[par]][k,], type="l", xlab="iteration", main=yname)
        }
      }
    }
  }
  return(g)
}
