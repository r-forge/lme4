fn <- "Coalition_lme4e_NM.RData"
source("Generic_lme4e_NM.R")

sessinfo <- sessionInfo()

data(coalition2,package="Zelig")

results <- lme4efit("Coalition",
                    coalition2,
                    duration ~ invest + fract + polar +
                       numst2 + crisis + (1 | country),
                    Gamma(link="log"))
save("sessinfo","results",file=fn)
