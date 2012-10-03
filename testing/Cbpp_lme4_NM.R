fn <- "Cbpp_lme4_NM.RData"
source("Generic_lme4.R")

sessinfo <- sessionInfo()
data(cbpp,package="lme4")
results <- lme4fit("Cbpp",
                   cbpp,
                   cbind(incidence,size-incidence)~period+(1|herd),
                   binomial)

save("results","sessinfo",file=fn)
