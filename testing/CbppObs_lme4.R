source("Generic_lme4.R")
fn <- "CbppObs_lme4.RData"
sessinfo <- sessionInfo()

data(cbpp)
cbpp$obs <- factor(seq(nrow(cbpp)))
results <- lme4fit("CbppObs",
                    cbpp,
                    cbind(incidence, size - incidence) ~ period + (1 | herd) + (1|obs),
                    binomial)

save("results","sessinfo",file=fn)

