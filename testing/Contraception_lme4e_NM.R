fn <- "Contraception_lme4e_NM.RData"
source("Generic_lme4e_NM.R")

sessinfo <- sessionInfo()

data(Contraception, package="mlmRev")
Contraception <- transform(Contraception,
                           ch=factor(ifelse(livch != 0, "Y", "N")))
results <- lme4efit("Contraception",
                    Contraception,
                    use ~ age + I(age^2) + ch + (1|district:urban),
                    binomial)
save("results",file=fn)
