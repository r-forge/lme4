library(lme4Eigen)

fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)

## 0.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
system.time( tpr <- profile(fm01ML) )
     
(confint(tpr) -> CIpr)
print(xyplot(tpr))
stopifnot(dim(CIpr) == c(3,2),
          all.equal(unname(CIpr[2,]),c(3.64362, 4.21446), tol=1e-6))


detach("package:lme4Eigen")
library(lme4a)
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)

## 0.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
system.time( tpr <- profile(fm01ML) )
     
(confint(tpr) -> CIpr)
print(xyplot(tpr))
stopifnot(dim(CIpr) == c(3,2),
          all.equal(unname(CIpr[2,]),c(3.64362, 4.21446), tol=1e-6))
