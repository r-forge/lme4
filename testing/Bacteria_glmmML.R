fn <- "Bacteria_glmmML.RData"
source("Generic_glmmML.R")

sessinfo <- sessionInfo()
data(bacteria,package="MASS")
bacteria$cluster <- bacteria$ID
if (is.factor(bacteria$y))
    bacteria$y <- as.numeric(bacteria$y)-1
form <- y~trt+I(week>2)
Family <- "binomial"

results <- glmmMLfit("Bacteria",bacteria,form,Family)

save("results","sessinfo",file=fn)


