fn <- "Contraception_NM1.RData"
## require(plyr)
source("xapply.R")
source("miscfuns.R")

library(lme4Eigen)
sessinfo <- sessionInfo()

data(Contraception, package="mlmRev")
Contraception <- within(Contraception,
                        ch <- factor(ifelse(livch != 0, "Y", "N")))

form <- use ~ age + I(age^2) + ch + (1|district:urban)

tolPwrssvec <- 10^seq(-11,-5,by=0.25)
nAGQvec <- c(2^(0:4),25)

results <- xapply(tolPwrssvec,nAGQvec,FUN=dfun_lme4Eigen_NM,
                  MoreArgs=list(problem="Contraception",
                  form=form,data=Contraception,family=binomial),
                  .progress="text",pbargs=list(style=3))
save("results",file=fn)
detach("package:lme4Eigen")

## testing: nAGQ=32?
## glmer(form=form,data=Contraception,family=binomial,nAGQ=25)
