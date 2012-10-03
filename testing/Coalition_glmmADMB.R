source("Generic_glmmADMB.R")
fn <- "Coalition_glmmADMB.RData"


sessinfo <- sessionInfo()

data(coalition2,package="Zelig")

results <- glmmADMBfit("Coalition",
                       coalition2,
                       duration ~ invest + fract + polar +
                       numst2 + crisis + (1 | country),
                    Family="Gamma",link="log",
                       extra.args="-ams 500000000",
                       verbose=TRUE)
save("sessinfo","results",file=fn)
