fn <- "Contraception_lme4e_bobyqa.RData"
source("Generic_lme4e.R")
library(optimx)
sessinfo <- sessionInfo()

data(Contraception, package="mlmRev")
Contraception <- transform(Contraception,
                           ch=factor(ifelse(livch != 0, "Y", "N")))
results <- lme4efit("Contraception",
                    method="bobyqa",
                    Contraception,
                    use ~ age + I(age^2) + ch + (1|district:urban),
                    binomial)
save("results","sessinfo",file=fn)
detach("package:optimx")
