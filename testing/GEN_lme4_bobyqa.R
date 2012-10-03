fn <- "NAME_lme4_bobyqa.RData"
source("Generic_lme4.R")
library(optimx)
sessinfo <- sessionInfo()

data(DATNAME, package="DATPKG")
results <- lme4fit("NAME",
                   method="bobyqa",
                   DATNAME,
                   FORM,
                   FAMILY)
save("results","sessinfo",file=fn)
detach("package:optimx")
