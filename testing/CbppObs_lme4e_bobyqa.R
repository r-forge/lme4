fn <- "CbppObs_lme4e_bobyqa.RData"
source("xapply.R")
source("miscfuns.R")

library(lme4Eigen)
library(optimx)
sessinfo <- sessionInfo()

data(cbpp)
cbpp$obs <- factor(seq(nrow(cbpp)))

form <- cbind(incidence, size - incidence) ~ period + (1 | herd) + (1|obs)

rhobegvec <- 2*10^(-8:0)
begendvec <- 10^(-3:-6)
tolPwrssvec <- 10^seq(-9,-5)

if (FALSE) { ## testing
  dfun_lme4Eigen_bobyqa(rhobegvec[1],begendvec[1],tolPwrssvec[1],form=form,data=cbpp,family=binomial,
                        problem="CbppObs")
}

results <- xapply(FUN=dfun_lme4Eigen_bobyqa,
                  rhobegvec,begendvec,tolPwrssvec,
                  MoreArgs=list(problem="CbppObs",
                    form=form,data=cbpp,family=binomial),
                  .progress="txt",pbargs=list(style=3))
save("results","sessinfo",file=fn)
