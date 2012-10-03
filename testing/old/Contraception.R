fn <- "Contraception1.RData"
source("xapply.R")
source("miscfuns.R")

library(lme4Eigen)
library(optimx)
sessinfo <- sessionInfo()

data(Contraception, package="mlmRev")
Contraception <- within(Contraception,
                        ch <- factor(ifelse(livch != 0, "Y", "N")))

form0 <- use ~ age + I(age^2) + ch + (1|district:urban)

d0 <- list(problem="Contraception",method="bobyqa",optimx=TRUE)
dbad <- c(d0,
          as.list(c(options=NA,time=NA,parameters=NA,
                    deviance=NA,KKT=NA,bad=NA,result=NA,
                    nAGQ=1)))
          
rhobegvec <- 2*10^(-8:0)
begendvec <- 10^(-3:-6)
tolPwrssvec <- 10^seq(-9,-5)

dfun <- function(rb,be,tolpwrss,...,debug=FALSE) {
  if (debug) cat(rb,be,tolpwrss,"\n")
  t0 <- system.time(cur.opt <- try(glmer(...,
                                         optimizer="bobyqa",
                                         tolPwrss=tolpwrss,
                                         control=list(rhobeg=rb,rhoend=rb*be)),
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
}
results <- xapply(FUN=dfun,rhobegvec,begendvec,tolPwrssvec,
                  MoreArgs=list(data=Contraception,form=form0,family=binomial),
                  .progress="txt",pbargs=list(style=3))
save("results",file=fn)
detach("package:lme4Eigen")
library(glmmADMB)

