source("Generic_lme4.R")
fn <- "Contraception_lme4.RData"
sessinfo <- sessionInfo()

data(Contraception, package="mlmRev")
Contraception <- transform(Contraception,
                           ch=factor(ifelse(livch != 0, "Y", "N")))

results <- lme4fit("Contraception",
                   Contraception,
                   use ~ age + I(age^2) + ch + (1|district:urban),
                   binomial)
save("results","sessinfo",file=fn)

