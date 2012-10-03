## utility functions for collecting bits of results
require(plyr)

## sapply with support for missing or broken values
Sapply <- function(x,fun,...) {
  r <- lapply(x,fun,...)
  r0 <- r[[1]]
  sapply(r,function(x) if (is.null(x) || is.null(x[[1]]) || is.na(x)) rep(NA,length(r0)) else x)
}

## get n'th element of each element in a list
getf <- function(x,n) {
  lapply(x,"[[",n)
}

pfun <- function(x) {
  paste(names(x),x,sep="=",collapse=",")
}

## should these data frames be assembled earlier, in the analysis functions?
getdat <- function(fn) {
  ## assemble data frame from results
  load(fn)
  n <- length(results)
  problem <- unlist(getf(results,"problem"))
  pkg <- unlist(getf(results,"pkg"))
  method <- unlist(getf(results,"method"))
  if (is.null(method)) method <- rep(NA,n)
  opts <- getf(results,"options")
  optlen <- max(c(1,sapply(opts,length)))
  opts2 <- laply(opts,
                 .drop=FALSE,
                 function(x) {
                   if (is.null(x) || all(is.na(x))) rep(NA,optlen) else unlist(x) })
  deviance <- Sapply(results,"[[","deviance")
  nAGQ <- Sapply(results,"[[","nAGQ")
  if (is.null(results[[1]]$KKT)) {
    KKT <- cbind(`1`=rep(NA,n),`2`=rep(NA,n))
  } else {
    KKT <- t(Sapply(results,"[[","KKT"))
  }
  time <- unlist(getf(results,"time"))
  pars <- do.call(rbind,
                  lapply(results,function(x) unlist(x$parameters)))
  dd <- data.frame(problem,pkg,method,opts2,deviance,KKT=KKT,pars,time,nAGQ)
  ## check that deviance is near minimum *within* calculation category
  id.vars <- c("opts2","tolPwrss","nAGQ","pkg","method")
  allid <- c(id.vars,c("rhobeg","rhoend"))
  dd2A <- ddply(dd,intersect(names(dd),id.vars),
                 function(x) {
                   scdev <- x$deviance-min(x$deviance,na.rm=TRUE)
                   devOK <- !is.na(scdev) & scdev<0.01
                   data.frame(x[,intersect(names(dd),allid)],devOK)
                 })
  dd2 <- merge(dd,dd2A)
  dd2$OKval <- with(dd2,ifelse(devOK & isTRUE(KKT.1) & isTRUE(KKT.2),
                             "good",
                             ifelse(devOK,"OK","bad")))
  if (nrow(dd2)!=n) stop("length mismatch")
  dd2
}


Rbind <- function(L) {
  ## merge a list of data frames
  allnames <- Reduce(union,lapply(L,names))
  L <- lapply(L,
              function(x) {
                extra.names <- setdiff(allnames,names(x))
                extra.cols <- matrix(NA,nrow=nrow(x),ncol=length(extra.names),
                                     dimnames=list(seq(nrow(x)),extra.names))
                data.frame(x,extra.cols)
              })
  do.call(rbind,L)
}
