source("xapply.R")
source("miscfuns.R")
library(plyr)
source("Generic_glmmML.R")

sessinfo <- sessionInfo()

fn <- "Contraception_glmmML.RData"
data(Contraception, package="mlmRev")
Contraception <- transform(Contraception,
                           ch=factor(ifelse(livch != 0, "Y", "N")),
                           use=as.numeric(use)-1,
                           cluster=interaction(district,urban))
form <- use ~ age + I(age^2) + ch
Family <- binomial

results <- glmmMLfit("Contraception",Contraception,form,binomial)

save("results","sessinfo",file=fn)
