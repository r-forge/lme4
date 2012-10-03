source("Generic_lme4.0.R")
fn <- "Cbpp_lme4.0.RData"
sessinfo <- sessionInfo()

data(cbpp,package="lme4")
results <- lme4.0fit("Cbpp",
                     cbpp,
                     cbind(incidence,size-incidence)~period+(1|herd),
                     binomial)

save("results","sessinfo",file=fn)
