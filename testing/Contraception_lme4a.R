source("Generic_lme4a.R")
fn <- "Contraception_lme4a.RData"
sessinfo <- sessionInfo()

data(Contraception, package="mlmRev")
Contraception <- transform(Contraception,
                           ch=factor(ifelse(livch != 0, "Y", "N")))

results <- lme4afit("Contraception",
                   Contraception,
                   use ~ age + I(age^2) + ch + (1|district:urban),
                   binomial)
save("results","sessinfo",file=fn)

