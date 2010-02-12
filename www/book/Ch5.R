###################################################
### chunk number 1: preliminaries
###################################################
options(width=74, show.signif.stars = FALSE,
        lattice.theme = function() canonical.theme("pdf", color = FALSE),
        str = strOptions(strict.width = "cut"))
library(splines)
library(lattice)
library(Matrix)
library(lme4a)


###################################################
### chunk number 2: fm8reduxshow eval=FALSE
###################################################
## fm8 <- lmer(Reaction ~ 1 + Days + (1 + Days|Subject), sleepstudy,
##             REML = 0, verbose = TRUE)


###################################################
### chunk number 3: fm8redux
###################################################
cat(paste(capture.output(fm8 <- lmer(Reaction ~ 1 + Days + (1 + Days|Subject),
                                     sleepstudy, REML = 0,
                                     verbose = TRUE))[1:4],
          collapse = "\n"), "\n...\n")


###################################################
### chunk number 4: attachrho
###################################################
attach(env(fm8))


###################################################
### chunk number 5: updateLambda
###################################################
str(Lambda)
str(Lind)
Lambda@x[] <- c(1,0,1)[Lind]
str(Lambda@x)
Ut <- crossprod(Lambda, Zt)


###################################################
### chunk number 6: updateL
###################################################
L <- update(L, Ut, mult = 1)


###################################################
### chunk number 7: RZXandRX
###################################################
RZX <- solve(L, solve(L, crossprod(Lambda, ZtX), sys = "P"), sys = "L")
RX <- chol(XtX - crossprod(RZX))


###################################################
### chunk number 8: cu
###################################################
cu <- solve(L, solve(L, crossprod(Lambda, Zty), sys = "P"), sys = "L")
cbeta <- solve(t(RX), Xty - crossprod(RZX, cu))


###################################################
### chunk number 9: stage2
###################################################
fixef <- as.vector(solve(RX, cbeta))
u <- solve(L, solve(L, cu - RZX %*% fixef, sys = "Lt"), sys = "Pt")


###################################################
### chunk number 10: wrapup
###################################################
mu <- gamma <- as.vector(crossprod(Ut, u) + X %*% fixef)
prss <- sum(c(y - mu, as.vector(u))^2)
ldL2 <- 2 * as.vector(determinant(L)$mod)
(deviance <- ldL2 + nobs * (1 + log(2 * pi * prss/nobs)))


###################################################
### chunk number 11: detachenv
###################################################
detach()


###################################################
### chunk number 12: cleanup
###################################################
rm(deviance, fixef, gamma, L, Lambda, ldL2, mu, prss, RX, RZX, u, Ut)


