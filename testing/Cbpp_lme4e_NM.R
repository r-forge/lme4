fn <- "Cbpp_lme4e_NM.RData"
source("Generic_lme4e_NM.R")

sessinfo <- sessionInfo()
data(cbpp,package="lme4Eigen")
results <- lme4efit("Cbpp",
                    cbpp,
                    cbind(incidence, size - incidence) ~ period + (1 | herd),
                    binomial)

save("results","sessinfo",file=fn)
