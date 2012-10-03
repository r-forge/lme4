source("Generic_lme4.0.R")
fn <- "NAME_lme4.0.RData"
sessinfo <- sessionInfo()

data(DATNAME,package="DATPKG")
results <- lme4.0fit("NAME",
                     DATNAME,
                     FORM,
                     FAMILY)

save("results","sessinfo",file=fn)
