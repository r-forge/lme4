fn <- "Cbpp_lme4_bobyqa.RData"
source("Generic_lme4.R")
library(optimx)
sessinfo <- sessionInfo()

data(cbpp, package="lme4")
results <- lme4fit("Cbpp",
                   method="bobyqa",
                   cbpp,
                   cbind(incidence,size-incidence)~period+(1|herd),
                   binomial)
save("results","sessinfo",file=fn)
detach("package:optimx")
