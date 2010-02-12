###################################################
### chunk number 1: preliminaries
###################################################
options(width=85, show.signif.stars = FALSE,
        lattice.theme = function() canonical.theme("pdf", color = FALSE),
        str = strOptions(strict.width = "cut"))
library(splines)
library(lattice)
library(Matrix)
library(lme4a)
data(Rail, package = "MEMSS")
fm1 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff)
fm1ML <- update(fm1, REML = FALSE)
if (file.exists("pr1.rda")) {
    load("pr1.rda")
} else {
    pr1 <- profile(fm1ML, delta = 0.2)
    save(pr1, file = "pr1.rda")
}


###################################################
### chunk number 2: lme4 eval=FALSE
###################################################
## library(lme4)


###################################################
### chunk number 3: strDyestuff
###################################################
str(Dyestuff)


###################################################
### chunk number 4: headDyestuff
###################################################
head(Dyestuff)


###################################################
### chunk number 5: summaryDyestuff
###################################################
summary(Dyestuff)


###################################################
### chunk number 6: Dyestuffdot
###################################################
set.seed(1234543)
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff,
              ylab = "Batch", jitter.y = TRUE, pch = 21,
              xlab = "Yield of dyestuff (grams of standard color)",
              type = c("p", "a")))


###################################################
### chunk number 7: strDye2
###################################################
str(Dyestuff2)
summary(Dyestuff2)


###################################################
### chunk number 8: Dyestuff2dot
###################################################
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff2,
              ylab = "Batch", jitter.y = TRUE, pch = 21,
              xlab = "Simulated response (dimensionless)",
              type = c("p", "a")))


###################################################
### chunk number 9: fm1
###################################################
fm1 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff)
print(fm1)


###################################################
### chunk number 10: fm1formula
###################################################
Yield ~ 1 + (1|Batch)


###################################################
### chunk number 11: fm1ML
###################################################
(fm1ML <- lmer(Yield ~ 1 + (1|Batch), Dyestuff, REML = FALSE))


###################################################
### chunk number 12: fm2
###################################################
(fm2 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff2))


###################################################
### chunk number 13: fm2ML
###################################################
(fm2ML <- update(fm2, REML = FALSE))


###################################################
### chunk number 14: fm2a
###################################################
summary(fm2a <- lm(Yield ~ 1, Dyestuff2))


###################################################
### chunk number 15: fm1MLverb
###################################################
fm1ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE, verbose = TRUE)


###################################################
### chunk number 16: fm1MLLambda
###################################################
env(fm1ML)$Lambda


###################################################
### chunk number 17: fm1Lambdaimage
###################################################
print(image(env(fm1)$Lambda, sub=NULL, xlab=NULL, ylab=NULL))


###################################################
### chunk number 18: fm1Ztimage
###################################################
print(image(env(fm1)$Zt, sub=NULL))


###################################################
### chunk number 19: pr1 eval=FALSE
###################################################
## pr1 <- profile(fm1ML)


###################################################
### chunk number 20: fm1prof
###################################################
print(xyplot(pr1, aspect = 1.3))


###################################################
### chunk number 21: fm1profshow
###################################################
xyplot(pr1, aspect = 1.3)


###################################################
### chunk number 22: confintfm195
###################################################
confint(pr1)


###################################################
### chunk number 23: confintfm199
###################################################
confint(pr1, level = 0.99)


###################################################
### chunk number 24: fm1absprof
###################################################
print(xyplot(pr1, absVal = TRUE, aspect = 0.7, strip=FALSE, strip.left = TRUE))


###################################################
### chunk number 25: sigmaprof
###################################################
zeta <- sqrt(qchisq(c(0.5,0.8,0.9,0.95,0.99), 1))
zeta <- c(-rev(zeta), 0, zeta)
spl <- attr(pr1, "forward")[[2]]
endpts <- predict(attr(pr1, "backward")[[2]], zeta)$y

fr <- data.frame(zeta = rep.int(zeta, 3),
                 endpts = c(endpts, exp(endpts), exp(2*endpts)),
                 pnm = gl(3, length(zeta)))
print(xyplot(zeta ~ endpts|pnm, fr, type = "h",
             scales = list(x = list(relation = "free")),
             xlab = NULL, ylab = expression(zeta), aspect = 1.3,
             strip = strip.custom(
             factor.levels = expression(log(sigma), sigma, sigma^2)),
             panel = function(...) {
                 panel.grid(h = -1, v = -1)
                 panel.abline(h=0)
                 panel.xyplot(...)
                 ll <- current.panel.limits()$xlim
                 lims <- switch(panel.number(), ll, log(ll), log(ll)/2)
                 pr <- predict(spl, seq(lims[1], lims[2], len = 101))
                 panel.lines(switch(panel.number(),
                                    pr$x,
                                    exp(pr$x),
                                    exp(pr$x * 2)), pr$y)
             }))


###################################################
### chunk number 26: sigma1prof
###################################################
zeta <- sqrt(qchisq(c(0.5,0.8,0.9,0.95,0.99), 1))
zeta <- c(-rev(zeta), 0, zeta)
spl <- attr(pr1, "forward")[[1]]
endpts <- predict(attr(pr1, "backward")[[1]], zeta)$y

fr <- data.frame(zeta = rep.int(zeta, 3),
                 endpts = c(log(endpts), endpts, endpts^2),
                 pnm = gl(3, length(zeta)))
## A mighty kludge here
fr[1,] <- c(NA, 1.5, 1)
fr[12,] <- c(NA, 0, 2)
print(xyplot(zeta ~ endpts|pnm, fr, type = "h",
             scales = list(x = list(relation = "free")),
             xlab = NULL, ylab = expression(zeta), aspect = 1.3,
             strip = strip.custom(
             factor.levels = expression(log(sigma[1]), sigma[1], sigma[1]^2)),
             panel = function(...) {
                 panel.grid(h = -1, v = -1)
                 panel.abline(h = 0)
                 panel.xyplot(...)
                 ll <- (current.panel.limits()$xlim)[2]
                 lims <- switch(panel.number(),
                                c(1.5, exp(ll)),
                                c(0, ll),
                                c(0, sqrt(ll)))
                 pr <- predict(spl, seq(lims[1], lims[2], len = 101))
                 panel.lines(switch(panel.number(),
                                    log(pr$x),
                                    pr$x,
                                    pr$x^2), pr$y)
             }))


###################################################
### chunk number 27: fm1profpair
###################################################
print(splom(pr1))


###################################################
### chunk number 28: fm1profpairshow eval=FALSE
###################################################
## splom(pr1)


###################################################
### chunk number 29: raneffm1
###################################################
ranef(fm1ML)


###################################################
### chunk number 30: strraneffm1
###################################################
str(ranef(fm1ML))


###################################################
### chunk number 31: fm1preddot
###################################################
print(dotplot(ranef(fm1ML, postVar=TRUE), strip = FALSE))


###################################################
### chunk number 32: fm1preddotshow eval=FALSE
###################################################
## dotplot(ranef(fm1ML, postVar = TRUE))


###################################################
### chunk number 33: fm1predqq
###################################################
print(qqmath(ranef(fm1ML, postVar=TRUE), strip = FALSE))


###################################################
### chunk number 34: fm1predqqshow eval=FALSE
###################################################
## qqmath(ranef(fm1ML, postVar=TRUE))


###################################################
### chunk number 35: Raildot
###################################################
print(dotplot(reorder(Rail, travel) ~ travel, Rail,
              ylab = "Rail", pch = 21,
              xlab = "Travel time for an ultrasonic wave (ms.)",
              type = c("p", "a")))


###################################################
### chunk number 36: attMEMSS eval=FALSE
###################################################
## library(MEMSS)


###################################################
### chunk number 37: loadRail eval=FALSE
###################################################
## data(Rail, package = "MEMSS")


###################################################
### chunk number 38: helpRail eval=FALSE
###################################################
## help(Rail, package = "MEMSS")


