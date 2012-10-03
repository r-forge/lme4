source("Generic_glmmADMB.R")
fn <- "NAME_glmmADMB.RData"

sessinfo <- sessionInfo()
data(DATNAME,package="DATPKG")

results <- glmmADMBfit("NAME",
                       DATNAME,
                       FORM,
                       FAMILY)
save("results","sessinfo",file=fn)


