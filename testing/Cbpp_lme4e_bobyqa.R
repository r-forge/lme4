fn <- "Cbpp_lme4e_bobyqa.RData"
source("xapply.R")
source("miscfuns.R")

source("Generic_lme4e.R")
library(optimx)
sessinfo <- sessionInfo()

data(cbpp, package="lme4Eigen")
results <- lme4efit("Cbpp",
                    method="bobyqa",
                    cbpp,
                    cbind(incidence, size - incidence) ~ period + (1 | herd),
                    binomial)
save("results","sessinfo",file=fn)
detach("package:optimx")
