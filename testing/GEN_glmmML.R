fn <- "NAME_glmmML.RData"
source("Generic_glmmML.R")

sessinfo <- sessionInfo()
data(DATNAME,package="DATPKG")
DATNAME$cluster <- DATNAME$CLUSTER
if (is.factor(DATNAME$RESPONSE))
    DATNAME$RESPONSE <- as.numeric(DATNAME$RESPONSE)-1
form <- FIXFORM
Family <- FAMILY

results <- glmmMLfit("NAME",DATNAME,form,Family)

save("results","sessinfo",file=fn)


