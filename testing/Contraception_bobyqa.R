fn <- "Contraception_bobyqa1.RData"
source("xapply.R")
source("miscfuns.R")

library(lme4Eigen)
library(optimx)
(sessinfo <- sessionInfo())

data(Contraception, package="mlmRev")
Contraception <- within(Contraception,
                        ch <- factor(ifelse(livch != 0, "Y", "N")))

form <- use ~ age + I(age^2) + ch + (1|district:urban)
        
rhobegvec <- 2*10^(-8:0)
begendvec <- 10^(-3:-6)
tolPwrssvec <- 10^seq(-9,-5)

if (FALSE) { ## testing
  dfun_lme4Eigen_bobyqa(rhobegvec[1],begendvec[1],tolPwrssvec[1],form=form,data=Contraception,family=binomial,
                        problem="Contraception")
}
results <- xapply(FUN=dfun_lme4Eigen_bobyqa,
                  rhobegvec,begendvec,tolPwrssvec,
                  MoreArgs=list(problem="Contraception",
                    form=form,data=Contraception,family=binomial),
                  .progress="txt",pbargs=list(style=3))
save("results",file=fn)
detach("package:lme4Eigen")

