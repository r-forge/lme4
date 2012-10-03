fn <- "Bacteria_lme4_NM.RData"
source("Generic_lme4.R")

sessinfo <- sessionInfo()
data(bacteria,package="MASS")
results <- lme4fit("Bacteria",
                   bacteria,
                   y~trt+I(week>2)+(1|ID),
                   "binomial")

save("results","sessinfo",file=fn)
