fn <- "Coalition_lme4.RData"
source("Generic_lme4.R")

sessinfo <- sessionInfo()

data(coalition2,package="Zelig")

results <- lme4fit("Coalition",
                    coalition2,
                    duration ~ invest + fract + polar +
                       numst2 + crisis + (1 | country),
                    Gamma(link="log"))
save("sessinfo","results",file=fn)
