
BFUN <- function(rhobeg,begend,tolPwrss,
                 form, Data, Family, npar, d0, dbad) {
  ff <- glmer(form, Data, Family, devFunOnly=2L)
  control <- list(rhobeg=rhobeg,rhoend=rhobeg*begend)
  ## tolPwrss is read from ff's environment
  tolPwrss0 <- get("tolPwrss",envir=environment(ff))
  assign("tolPwrss",tolPwrss,envir=environment(ff))
  t0 <- system.time(cur.opt <- try(optimx(val, ff,
                                          lower=c(0, rep(-Inf, npar)),
                                          control=control,
                                          method=method)))
  assign("tolPwrss",tolPwrss0,envir=environment(ff))
  if (inherits(cur.opt,"try-error")) {
    d <- dbad
    d$result <- attr(cur.opt,"condition")$message
  } else {
    d <- d0
    d <- c(d,list(options=c(control,list(tolPwrss=tolPwrss)),
                  time=sum(t0[1:2]), ## user +system
                  parameters=cur.opt$par$par,           ## ??? why nested
                  deviance=c(cur.opt$fvalues$fvalues),  ## ??? why nested
                  KKT=c(cur.opt$KKT1,cur.opt$KKT2),
                  bad=NA,
                  result=NA_character_))
  }
}

xapply <- function(FUN,...,FLATTEN=TRUE,MoreArgs=NULL) {
  ## add progress bar??
  L <- list(...)
  inds <- do.call(expand.grid,lapply(L,seq_along)) ## Marek's suggestion
  retlist <- vector("list",nrow(inds))
  for (i in 1:nrow(inds)) {
    arglist <- mapply(function(x,j) x[[j]],L,as.list(inds[i,]),SIMPLIFY=FALSE)
    if (FLATTEN) {
      retlist[[i]] <- do.call(FUN,c(arglist,MoreArgs))
    }
  }
  retlist
}



do_bobyqa <- function(method="bobyqa",Data=Contraception,
                      form=use ~ age + I(age^2) + ch + (1|district:urban),
                      Family=binomial,
                      parameters=list(rhobeg=2*10^(-8:0),
                        begend=10^(-3:-6),
                        tolPwrss=10^seq(-9,-5))) {
  dname <- deparse(substitute(Data))
  sessinfo <- sessionInfo()
  f0 <- glmer(form, Data, Family)
  ## FIXME: maybe overkill? All we need at this point is
  ##   ncol(model.matrix(~form,Data)) [but RE spec. will get in the way]
  nbeta <- length(fixef(f0)) 
  ## extract the objective function for the second stage optimization
  ff <- glmer(form, Data, Family, devFunOnly=2L)
  val <- c(get("pp", environment(ff))$theta, get("pp", environment(ff))$beta(0))
  d0 <- list(problem=dname,method=method,optimx=TRUE)
  dbad <- c(d0,as.list(c(options=NA,time=NA,parameters=NA,
                         deviance=NA,KKT=NA,bad=NA,result=NA)))
  do.call(xapply,
          parameters,
          FUN=BFUN,MoreArgs=list(form=form,
                     Data=Data,Family=Family,npar=nbeta,
                     d0=d0,dbad=dbad))
}

library(lme4Eigen)
library(optimx)

data(Contraception, package="mlmRev")
Contraception <- within(Contraception,
                        ch <- factor(ifelse(livch != 0, "Y", "N")))

do_bobyqa()
do_bobyqa <- function(method="bobyqa",Data=Contraception,
                      form=use ~ age + I(age^2) + ch + (1|district:urban),
                      Family=binomial,
                      parameters=list(rhobeg=2*10^(-8:0),
                        begend=10^(-3:-6),
                        tolPwrss=10^seq(-9,-5)))
