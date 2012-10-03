fn <- "Contraception_bobyqa1.RData"
source("../xapply.R")
source("../miscfuns.R")

library(lme4Eigen)
library(optimx)
sessinfo <- sessionInfo()

data(Contraception, package="mlmRev")
Contraception <- within(Contraception,
                        ch <- factor(ifelse(livch != 0, "Y", "N")))

form <- use ~ age + I(age^2) + ch + (1|district:urban)


results <- list()
library(lme4)
t1 <- system.time(g1 <- glmer(form,Contraception,binomial))
results <- c(list(problem="Contraception",
                  pkg="lme4",method="nlminb",optimx=FALSE)
