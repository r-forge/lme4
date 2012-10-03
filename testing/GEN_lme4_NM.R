fn <- "NAME_lme4_NM.RData"
source("Generic_lme4.R")

sessinfo <- sessionInfo()
data(DATNAME,package="DATPKG")
results <- lme4fit("NAME",
                   DATNAME,
                   FORM,
                   FAMILY)

save("results","sessinfo",file=fn)
