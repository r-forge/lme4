fn <- "Bacteria_lme4_bobyqa.RData"
source("Generic_lme4.R")
library(optimx)
sessinfo <- sessionInfo()

data(bacteria, package="MASS")
results <- lme4fit("Bacteria",
                   method="bobyqa",
                   bacteria,
                   y~trt+I(week>2)+(1|ID),
                   "binomial")
save("results","sessinfo",file=fn)
detach("package:optimx")
