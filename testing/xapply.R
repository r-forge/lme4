.simpleCap <- function(x) {
  paste(toupper(substring(x, 1,1)), substring(x, 2),
        sep="", collapse=" ")
}

xapply <- function(FUN,...,FLATTEN=TRUE,MoreArgs=NULL,
                   namelist,.progress="none",pbargs=list()) {
  L <- list(...)
  if (.progress=="text") .progress <- "txt" ## compatibility with plyr
  inds <- do.call(expand.grid,lapply(L,seq_along)) ## Marek's suggestion
  maxind <- nrow(inds)
  if (.progress!="none") {  ## progress bar
    pbfun <- get(paste(.progress,"ProgressBar",sep=""))
    setpbfun <- get(paste("set",.simpleCap(.progress),"ProgressBar",sep=""))
    pb <- do.call(pbfun,pbargs)
    FUN0 <- FUN
    FUN <- function(...) {
      setpbfun(pb,i/maxind)
      FUN0(...)
    }
  }
  retlist <- vector("list",nrow(inds))
  for (i in 1:nrow(inds)) {
    arglist <- mapply(function(x,j) x[[j]],L,as.list(inds[i,]),SIMPLIFY=FALSE)
    if (FLATTEN) {
      retlist[[i]] <- do.call(FUN,c(arglist,MoreArgs))
    }
  }
  if (!missing(namelist)) attr(retlist,"grid") <- do.call(expand.grid,namelist)
  if (.progress!="none") close(pb)
  retlist
}

if (FALSE) {
  ## test example
  L1 <- list(data.frame(x=1:10,y=1:10),
             data.frame(x=runif(10),y=runif(10)),
             data.frame(x=rnorm(10),y=rnorm(10)))
  L2 <- list(y~1,y~x,y~poly(x,2))          
  z <- xapply(lm,L2,L1)
  require(tcltk)
  z <- xapply(lm,L2,L1,.progress="tk")
}
