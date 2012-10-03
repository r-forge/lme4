fn <- "Cbpp_glmmML.RData"
source("Generic_glmmML.R")

sessinfo <- sessionInfo()
data(cbpp,package="lme4")
cbpp$cluster <- cbpp$herd
if (is.factor(cbpp$incidence))
    cbpp$incidence <- as.numeric(cbpp$incidence)-1
form <- cbind(incidence,size-incidence)~period
Family <- binomial

results <- glmmMLfit("Cbpp",cbpp,form,Family)

save("results","sessinfo",file=fn)


