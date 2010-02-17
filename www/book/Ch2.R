###################################################
### chunk number 1: preliminaries
###################################################
options(lattice.theme = function() canonical.theme("pdf", color = FALSE),
        str = strOptions(strict.width = "cut"))
library(splines)
library(lattice)
library(Matrix)
library(lme4a)
fm2ML <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin, REML=0)
fm3 <- lmer(strength ~ 1 + (1|sample) + (1|batch), Pastes, REML=0)
fm3a <- lmer(strength ~ 1 + (1|sample), Pastes, REML=0)
if (file.exists("pr2.rda")) {
    load("pr2.rda")
} else {
    pr2 <- profile(fm2ML)
    save(pr2, file = "pr2.rda")
}
if (file.exists("pr3.rda")) {load("pr3.rda")
} else {
    pr3 <- profile(fm3)
    save(pr3, file = "pr3.rda")
}
if (file.exists("pr3a.rda")) {load("pr3a.rda")
} else {
    pr3a <- profile(fm3a)
    save(pr3a, file = "pr3a.rda")
}
if (file.exists("fm4.rda")) {load("fm4.rda")
} else {
    fm4 <- lmer(y ~ 1 + (1|s) + (1|d) + (1|dept:service), InstEval, REML=0)
    save(fm4, file = "fm4.rda")
}
if (file.exists("rr4.rda")) {load("rr4.rda")
} else {                             
    rr4 <- ranef(fm4, postVar = TRUE, whichel = "dept:service")
    save(rr4, file = "rr4.rda")
}
if (file.exists("fm4a.rda")) { load("fm4a.rda")
} else {
    fm4a <- lmer(y ~ 1 + (1|s) + (1|d), InstEval, REML=0)
    save(fm4a, file = "fm4a.rda")
}
if (!file.exists("figs/Multiple-fm4Limage.png")) {
    ttt <- env(fm4)$L
    ttt@x[] <- 1
    trellis.device(png, color = FALSE, file = "figs/Multiple-fm4Limage.png",
                   height=6, width=8, units='in', res=600)
    print(image(ttt, sub=NULL, xlab=NULL, ylab=NULL, colorkey=FALSE))
    rm(ttt)
    dev.off()
}
options(width=74, show.signif.stars = FALSE,
        str = strOptions(strict.width = "cut"))


###################################################
### chunk number 2: strPen
###################################################
str(Penicillin)


###################################################
### chunk number 3: sumPen
###################################################
summary(Penicillin)


###################################################
### chunk number 4: Penicillindot
###################################################
print(dotplot(reorder(plate, diameter) ~ diameter, Penicillin, groups = sample,
              ylab = "Plate", xlab = "Diameter of growth inhibition zone (mm)",
              type = c("p", "a"), auto.key = list(columns = 3, lines = TRUE)))


###################################################
### chunk number 5: xtabsPenicillin
###################################################
xtabs(~ sample + plate, Penicillin)


###################################################
### chunk number 6: fm2
###################################################
(fm2 <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin))


###################################################
### chunk number 7: fm2ranef
###################################################
qrr2 <- dotplot(ranef(fm2, postVar = TRUE), strip = FALSE)
print(qrr2[[1]], pos = c(0,0,1,0.75), more = TRUE)
print(qrr2[[2]], pos = c(0,0.65,1,1))


###################################################
### chunk number 8: fm2Ztimage
###################################################
print(image(env(fm2)$Zt, sub=NULL,xlab=NULL,ylab=NULL))


###################################################
### chunk number 9: fm2LambdaLimage
###################################################
print(image(env(fm2)$Lambda, sub=NULL,xlab=expression(Lambda),ylab=NULL),
      split=c(1,1,3,1), more=TRUE)
print(image(tcrossprod(env(fm2)$Zt),sub=NULL,xlab="Z'Z",ylab=NULL),
      split=c(2,1,3,1), more=TRUE)
print(image(env(fm2)$L, sub=NULL,xlab="L",ylab=NULL,colorkey=FALSE),
      split=c(3,1,3,1))


###################################################
### chunk number 10: fm2verb
###################################################
env(fm2)$theta


###################################################
### chunk number 11: fm2prplot
###################################################
print(xyplot(pr2, absVal=0, aspect=1.3, layout=c(4,1)))


###################################################
### chunk number 12: fm2confint
###################################################
confint(pr2)


###################################################
### chunk number 13: fm2confintrel
###################################################
confint(pr2)[1:2,]/c(0.8455722, 1.770648)


###################################################
### chunk number 14: fm2prpairs
###################################################
print(splom(pr2))


###################################################
### chunk number 15: fm2lprpairs
###################################################
print(splom(log(pr2)))


###################################################
### chunk number 16: strsumPastes
###################################################
str(Pastes)
summary(Pastes)


###################################################
### chunk number 17: imagextabsPastes
###################################################
print(image(xtabs(~ batch + sample, Pastes, sparse = TRUE),
            sub = NULL, xlab = "sample", ylab = "batch"))


###################################################
### chunk number 18: xtabsPastes
###################################################
xtabs(~ batch + sample, Pastes, drop = TRUE, sparse = TRUE)


###################################################
### chunk number 19: Pastesplot
###################################################
pp <- Pastes
pp <- within(pp, bb <- reorder(batch, strength))
pp <- within(pp, ss <- reorder(reorder(sample, strength),
          as.numeric(batch)))
print(dotplot(ss ~ strength | bb, pp, pch = 21,
              strip = FALSE, strip.left = TRUE, layout = c(1, 10),
              scales = list(y = list(relation = "free")),
              ylab = "Sample within batch", type = c("p", "a"),
              xlab = "Paste strength", jitter.y = TRUE))


###################################################
### chunk number 20: Pastesxtab
###################################################
xtabs(~ cask + batch, Pastes)


###################################################
### chunk number 21: nested1
###################################################
Pastes$sample <- with(Pastes, factor(batch:cask))


###################################################
### chunk number 22: nested2
###################################################
Pastes <- within(Pastes, sample <- factor(batch:cask))


###################################################
### chunk number 23: fm3
###################################################
(fm3 <- lmer(strength ~ 1 + (1|sample) + (1|batch), Pastes, REML=0))


###################################################
### chunk number 24: fm3LambdaLimage
###################################################
print(image(env(fm3)$Lambda, sub=NULL,xlab=expression(Lambda),ylab=NULL),
      split=c(1,1,3,1), more=TRUE)
print(image(tcrossprod(env(fm3)$Zt),sub=NULL,xlab="Z'Z",ylab=NULL),
      split=c(2,1,3,1), more=TRUE)
print(image(env(fm3)$L, sub=NULL,xlab="L",ylab=NULL,colorkey=FALSE),
      split=c(3,1,3,1))


###################################################
### chunk number 25: stddev
###################################################
f3 <- "%.3f"
f4 <- "%.4f"
vc3 <- VarCorr(fm3)
stddev <- unname(c(sapply(vc3, function(el) attr(el, "stddev")),
                   attr(vc3, "sc")))


###################################################
### chunk number 26: fm3ranef
###################################################
qrr3 <- dotplot(ranef(fm3, postVar = TRUE), strip = FALSE)
print(qrr3[[1]], pos = c(0,0,1,0.75), more = TRUE)
print(qrr3[[2]], pos = c(0,0.65,1,1))


###################################################
### chunk number 27: fm3prplot
###################################################
print(xyplot(pr3, absVal=0, aspect=1.3, layout=c(4,1)))


###################################################
### chunk number 28: fm3a
###################################################
(fm3a <- lmer(strength ~ 1 + (1|sample), Pastes, REML=0)) 


###################################################
### chunk number 29: anovafm3
###################################################
anova(fm3a, fm3)


###################################################
### chunk number 30: fm3aprplot
###################################################
print(xyplot(pr3a, aspect=1.3, layout=c(3,1)))


###################################################
### chunk number 31: confintfm3
###################################################
confint(pr3)
confint(pr3a)


###################################################
### chunk number 32: fm3aprpairs
###################################################
print(splom(pr3a))


###################################################
### chunk number 33: strInstEval
###################################################
## the default  strict.width is "no" for back-compatibility
## instead, we could also *globally* set
##    strOptions(strict.width = "cut")
##str(InstEval, strict.width = "cut")
str(InstEval)


###################################################
### chunk number 34: xtabsy
###################################################
xtabs(~ y, InstEval)


###################################################
### chunk number 35: fm4show eval=FALSE
###################################################
## (fm4 <- lmer(y ~ 1 + (1|s) + (1|d)+(1|dept:service), InstEval, REML=0))


###################################################
### chunk number 36: fm4
###################################################
fm4


###################################################
### chunk number 37: fm4ranef
###################################################
print(dotplot(rr4, strip = FALSE))


###################################################
### chunk number 38: LRTshow eval=FALSE
###################################################
## fm4a <- lmer(y ~ 1 + (1|s) + (1|d), InstEval, REML=0)
## anova(fm4a,fm4)


###################################################
### chunk number 39: LRTfm4
###################################################
anova(fm4a,fm4)


###################################################
### chunk number 40: fm4Lstats
###################################################
object.size(env(fm4)$L)
unclass(round(object.size(env(fm4)$L)/2^20, 3))  # size in megabytes


###################################################
### chunk number 41: fm4densestats
###################################################
(8 * (4128 * 4129)/2)/2^20   # size in megabytes


###################################################
### chunk number 42: fm4densestats
###################################################
(8 * 4128^2)/2^20   # size in megabytes


###################################################
### chunk number 43: nnzero
###################################################
nnzero(as(env(fm4)$L, "sparseMatrix"))


###################################################
### chunk number 44: badformula
###################################################
strength ~ 1 + (1|cask) + (1|batch)


###################################################
### chunk number 45: slashform
###################################################
strength ~ 1 + (1|batch/cask)


###################################################
### chunk number 46: expandedform
###################################################
strength ~ 1 + (1|batch) + (1|batch:cask)


###################################################
### chunk number 47: attMEMSS eval=FALSE
###################################################
## library(MEMSS)


###################################################
### chunk number 48: loadergo eval=FALSE
###################################################
## data(ergoStool, package = "MEMSS")


###################################################
### chunk number 49: xtabsergo eval=FALSE
###################################################
## xtabs(~ Type + Subject, ergoStool)


###################################################
### chunk number 50: plotexam eval=FALSE
###################################################
## dotplot(ranef(fm, which = "Type", postVar = TRUE), aspect = 0.2,
##         strip = FALSE)


###################################################
### chunk number 51: xtabsPairsStorage eval=FALSE
###################################################
## xtabs(~ Pair + Storage, Meat)


