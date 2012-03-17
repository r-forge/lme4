library(lme4.0)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
stopifnot(all.equal(attr(VarCorr(fm1),"sc"),sigma(fm1)))
