source("Generic_glmmADMB.R")
fn <- "Contraception_glmmADMB.RData"

sessinfo <- sessionInfo()

data(Contraception, package="mlmRev")
Contraception <- within(Contraception,
                        ch <- factor(ifelse(livch != 0, "Y", "N")))

results <- glmmADMBfit("Contraception",
                       Contraception,
                       Form=use ~ age + I(age^2) + ch + (1|district:urban),
                       binomial)
save("results","sessinfo",file=fn)


