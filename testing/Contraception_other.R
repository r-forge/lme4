fn <- "Contraception_lme4.RData"
## require(plyr)
source("xapply.R")
source("miscfuns.R")

library(lme4)
sessinfo <- sessionInfo()

data(Contraception, package="mlmRev")
Contraception <- within(Contraception,
                        ch <- factor(ifelse(livch != 0, "Y", "N")))

form <- use ~ age + I(age^2) + ch + (1|district:urban)

nAGQvec <- 2^(0:5)

results <- xapply(nAGQvec,FUN=dfun_lme4other,
                  MoreArgs=list(problem="Contraception",pkg="lme4",
                    form=form,data=Contraception,family=binomial),
                  .progress="text",pbargs=list(style=3))
detach("package:lme4",unload=TRUE)
save("results",file=fn)
library(lme4a)

## lme4a does not do AGQ ...
results <- c(results,list(dfun_lme4other(nAGQ=1,
                                problem="Contraception",pkg="lme4a",
                                form=form,data=Contraception,family=binomial)))
detach("package:lme4a",unload=TRUE)
## Warning messages:
## 1: In FUN(X[[2L]], ...) :
##   Created a package name, ‘2012-01-21 16:05:59’, when none found

save("results",file=fn)


