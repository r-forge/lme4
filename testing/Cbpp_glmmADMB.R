source("Generic_glmmADMB.R")
fn <- "Cbpp_glmmADMB.RData"

sessinfo <- sessionInfo()
data(cbpp,package="lme4")

results <- glmmADMBfit("Cbpp",
                       cbpp,
                       cbind(incidence,size-incidence)~period+(1|herd),
                       binomial)
save("results","sessinfo",file=fn)


