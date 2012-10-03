source("Generic_glmmADMB.R")
fn <- "Bacteria_glmmADMB.RData"

sessinfo <- sessionInfo()
data(bacteria,package="MASS")

results <- glmmADMBfit("Bacteria",
                       bacteria,
                       y~trt+I(week>2)+(1|ID),
                       "binomial")
save("results","sessinfo",file=fn)


