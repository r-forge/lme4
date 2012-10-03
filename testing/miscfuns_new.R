## FIXME: merge with coeftab?  Am I reinventing it?
allcoef0 <- function(x,w=c("fixef","VarCorr")) {  ## "ranef" also OK
  res <- list()
  if ("fixef" %in% w) {
    if (inherits(x,"mer") || inherits(x,"merMod") || inherits(x,"glmmadmb")) {
      res <- c(res,list(fixef=fixef(x)))
    }
  }
  if ("ranef" %in% w) {
    res <- c(res,list(ranef=unlist(ranef(x))))
  }
  if ("VarCorr" %in% w) {
    ## FIXME: test for complex V-C spec
    vc <- VarCorr(x)
    if (inherits(x,"mer") || inherits(x,"merMod")) {
      vpars <- sapply(vc,
                      function(x) {
                        vv <- attr(x,"stddev")
                        if (nrow(x)>1) {
                          corr <- attr(x,"correlation")
                          vv <- c(vv,corr[lower.tri(corr)])
                        }
                        vv
                      })
      sc <- attr(vc,"sc")
      if (!is.na(sc) && sc!=1) {
        vpars <- c(vpars,sc)
      }
    } else {
      ## FIXME
      vpars <- sqrt(unlist(vc))
    }
    res <- c(res,list(VarCorr=vpars))
  }
  res
}

require(coefplot2)
allcoef <- function(x) {
  tmpf <- function(x) {
    r <- x[,"Estimate"]
    names(r) <- rownames(x)
    r
  }
  list(fixef=tmpf(coeftab(x,ptype="fixef")),
       VarCorr=tmpf(coeftab(x,ptype="vcov")))
}

dfun_lme4Eigen_bobyqa <- function(rb,be,tolpwrss,nAGQ=1,...,problem,debug=FALSE) {
  if (debug) cat(rb,be,tolpwrss,"\n")
  control <- list(rhobeg=rb,rhoend=rb*be)
  d0 <- list(problem=problem,pkg="lme4Eigen",method="bobyqa",optimx=TRUE,nAGQ=nAGQ,
             options=c(control,list(tolPwrss=tolpwrss)))
  dbad <- c(d0,
            as.list(c(time=NA,parameters=NA,
                      deviance=NA,KKT=NA,bad=NA,result=NA)))
  t0 <- system.time(cur.opt <- try(glmer(...,nAGQ=nAGQ,
                                         optimizer="bobyqa",
                                         tolPwrss=tolpwrss,
                                         control=control),
                                   silent=TRUE))
  if (inherits(cur.opt,"try-error")) {
    d <- dbad
    d$result <- attr(cur.opt,"condition")$message
  } else {
    d <- d0
    d <- c(d,list(options=c(control,list(tolPwrss=tolpwrss)),
                  time=sum(t0[1:2]), ## user +system
                  parameters=allcoef(cur.opt),
                  deviance=deviance(cur.opt),
                  KKT=c(NA,NA),
                  ## FIXME: extract KKT from lme4Eigen?
                  bad=NA,
                  result=NA_character_))
  }
  d
}

dfun_lme4Eigen_NM <- function(tolpwrss,nAGQ,...,problem,debug=FALSE) {
  d0 <- list(problem=problem,pkg="lme4Eigen",method="NelderMead",optimx=FALSE,nAGQ=nAGQ,
             options=c(list(tolPwrss=tolpwrss)))
  dbad <- c(d0,
            as.list(c(time=NA,parameters=NA,
                      deviance=NA,KKT=NA,bad=NA,result=NA)))
  t0 <- system.time(cur.opt <- try(glmer(...,
                                         optimizer="NelderMead",
                                         nAGQ=nAGQ,
                                         tolPwrss=tolpwrss),
                                   silent=TRUE))
  if (inherits(cur.opt,"try-error")) {
    d <- dbad
    d$result <- attr(cur.opt,"condition")$message
  } else {
    d <- d0
    d <- c(d,list(options=list(tolPwrss=tolpwrss),
                  time=sum(t0[1:2]), ## user +system
                  parameters=allcoef(cur.opt),
                  deviance=deviance(cur.opt),
                  KKT=c(NA,NA),
                  bad=NA,
                  nAGQ=nAGQ,
                  result=NA_character_))
  }
  d
}

dfun_lme4other <- function(nAGQ,...,problem,pkg,debug=FALSE) {
  if (debug) cat(problem,pkg,nAGQ,"\n")
  d0 <- list(problem=problem,pkg=pkg,
             method="nlminb",optimx=FALSE,nAGQ=nAGQ,options=NULL)
  t0 <- system.time(cur.opt <- try(glmer(...,nAGQ=nAGQ),
                                   silent=TRUE))
  if (inherits(cur.opt,"try-error")) {
    d <- dbad
    d$result <- attr(cur.opt,"condition")$message
  } else {
    d <- d0
    d <- c(d,list(options=NULL,
                  time=sum(t0[1:2]), ## user +system
                  parameters=allcoef(cur.opt),
                  deviance=deviance(cur.opt),
                  KKT=c(NA,NA),
                  bad=NA,
                  nAGQ=nAGQ,
                  result=NA_character_))
  }
  d
}
