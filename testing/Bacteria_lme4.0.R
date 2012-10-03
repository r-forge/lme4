source("Generic_lme4.0.R")
fn <- "Bacteria_lme4.0.RData"
sessinfo <- sessionInfo()

data(bacteria,package="MASS")
results <- lme4.0fit("Bacteria",
                     bacteria,
                     y~trt+I(week>2)+(1|ID),
                     "binomial")

save("results","sessinfo",file=fn)
