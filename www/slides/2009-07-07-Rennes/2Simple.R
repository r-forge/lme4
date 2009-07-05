###################################################
### chunk number 1: preliminaries
###################################################
options(width=69,show.signif.stars=FALSE)
library(lattice)
library(Matrix)
library(lme4)
#lattice.options(default.theme = function() standard.theme())
lattice.options(default.theme = function() standard.theme(color=FALSE))
if (exists("classroom.rda")) {
    load("classroom.rda")
} else {
    classroom <- within(read.csv("http://www-personal.umich.edu/~bwest/classroom.csv"),
                    {
                        classid <- factor(classid)
                        schoolid <- factor(schoolid)
                        sex <- factor(sex, labels = c("M","F"))
                        minority <- factor(minority, labels = c("N", "Y"))
                    })
}


###################################################
### chunk number 2: installlme4 eval=FALSE
###################################################
## install.packages("lme4")


###################################################
### chunk number 3: require eval=FALSE
###################################################
## require(lme4)


###################################################
### chunk number 4: attach eval=FALSE
###################################################
## library(lme4)


###################################################
### chunk number 5: Dyestuffstr
###################################################
str(Dyestuff)
summary(Dyestuff)


###################################################
### chunk number 6: Dyestuffplot
###################################################
set.seed(1234543)
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff,
              ylab = "Batch", jitter.y = TRUE, pch = 21, aspect = 0.32,
              xlab = "Yield of dyestuff (grams of standard color)",
              type = c("p", "a")))


###################################################
### chunk number 7: fm1
###################################################
fm1 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff)
print(fm1)


###################################################
### chunk number 8: op
###################################################
op <- options(digits=5)


###################################################
### chunk number 9: extractors
###################################################
fixef(fm1)
ranef(fm1, drop = TRUE)
fitted(fm1)


###################################################
### chunk number 10: unop
###################################################
options(op)


###################################################
### chunk number 11: strfm1
###################################################
str(model.matrix(fm1))
fm1@Zt


###################################################
### chunk number 12: fm1ranef
###################################################
print(dotplot(ranef(fm1, postVar = TRUE), strip = FALSE)[[1]])


###################################################
### chunk number 13: update
###################################################
(fm1M <- update(fm1, REML = FALSE))


###################################################
### chunk number 14: fm1refit
###################################################
invisible(update(fm1, verbose = TRUE))


###################################################
### chunk number 15: Dyestuff2
###################################################
str(Dyestuff2)


###################################################
### chunk number 16: Dyestuff2plot
###################################################
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff2,
              ylab = "Batch", jitter.y = TRUE, pch = 21, aspect = 0.42,
              xlab = "Simulated response (dimensionless)",
              type = c("p", "a")))


###################################################
### chunk number 17: fm1A
###################################################
(fm1A <- lmer(Yield ~ 1 + (1|Batch), Dyestuff2, verbose = TRUE))


###################################################
### chunk number 18: lm1
###################################################
summary(lm1 <- lm(Yield ~ 1, Dyestuff2))
logLik(lm1)


###################################################
### chunk number 19: fm1call
###################################################
fm1@call


###################################################
### chunk number 20: Penicillinstr
###################################################
str(Penicillin)
xtabs(~ sample + plate, Penicillin)


###################################################
### chunk number 21: PenicillinPlot
###################################################
print(dotplot(reorder(plate, diameter) ~ diameter, Penicillin, groups = sample,
              ylab = "Plate", xlab = "Diameter of growth inhibition zone (mm)",
              type = c("p", "a"), auto.key = list(columns = 6, lines = TRUE)))


###################################################
### chunk number 22: fm2
###################################################
(fm2 <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin))


###################################################
### chunk number 23: 
###################################################
op <- options(digits = 5)


###################################################
### chunk number 24: fixef2
###################################################
fixef(fm2)
ranef(fm2, drop = TRUE)


###################################################
### chunk number 25: 
###################################################
options(op)


###################################################
### chunk number 26: fm2ranef
###################################################
qrr2 <- dotplot(ranef(fm2, postVar = TRUE), strip = FALSE)
print(qrr2[[1]], pos = c(0,0,1,0.75), more = TRUE)
print(qrr2[[2]], pos = c(0,0.65,1,1))


###################################################
### chunk number 27: fm2Z
###################################################
print(image(fm2@Zt, xlab = NULL, ylab = NULL, sub = "Z'"))


###################################################
### chunk number 28: fm2formula
###################################################
fm2@call[["formula"]]


###################################################
### chunk number 29: Pastesstr
###################################################
str(Pastes)
xtabs(~ batch + sample, Pastes, sparse = TRUE)


###################################################
### chunk number 30: samplegen eval=FALSE
###################################################
## Pastes <- within(Pastes, sample <- factor(batch:cask))


###################################################
### chunk number 31: Pastesplot
###################################################
Pastes <- within(Pastes, bb <- reorder(batch, strength))
Pastes <- within(Pastes, ss <- reorder(reorder(sample, strength),
          as.numeric(batch)))
print(dotplot(ss ~ strength | bb, Pastes,
              strip = FALSE, strip.left = TRUE, layout = c(1, 10),
              scales = list(y = list(relation = "free")),
              ylab = "Sample within batch", type = c("p", "a"),
              xlab = "Paste strength", jitter.y = TRUE))


###################################################
### chunk number 32: fm3
###################################################
(fm3 <- lmer(strength ~ 1 + (1|batch) + (1|sample), Pastes))


###################################################
### chunk number 33: fm3ranef
###################################################
qrr3 <- dotplot(ranef(fm3, postVar = TRUE), strip = FALSE)
print(qrr3[[1]], pos = c(0,0,1,0.75), more = TRUE)
print(qrr3[[2]], pos = c(0,0.65,1,1))


###################################################
### chunk number 34: Z3fig
###################################################
print(image(fm3@Zt, xlab = NULL, ylab = NULL, sub = NULL))


###################################################
### chunk number 35: fm3LRT
###################################################
fm3M <- update(fm3, REML = FALSE)
fm4M <- lmer(strength ~ 1 + (1|sample), Pastes, REML = FALSE)
anova(fm4M, fm3M)


###################################################
### chunk number 36: fm4
###################################################
(fm4 <- update(fm4M, REML = TRUE))


###################################################
### chunk number 37: xtabsclass
###################################################
xtabs( ~ xtabs(~ classid, classroom))


###################################################
### chunk number 38: xtabsclass
###################################################
table(xtabs(~ schoolid,
            unique(subset(classroom, select = c(classid, schoolid)))))


###################################################
### chunk number 39: Schoolsplot
###################################################
refactor <- function(x) if(is.factor(x)) factor(x) else x
sch12 <- do.call(data.frame,
                 lapply(subset(classroom,
                               schoolid %in% c(12,15, 17, 33,46, 57,
                                               68, 70, 71, 76, 85, 99)),
                        refactor))
sch12 <- within(sch12, ss <- reorder(schoolid, mathgain))
sch12 <- within(sch12, cc <- reorder(reorder(classid, mathgain),
          as.numeric(schoolid)))
print(dotplot(cc ~ mathgain | ss , sch12, 
              strip = FALSE, strip.left = TRUE, layout = c(1, 12),
              scales = list(y = list(relation = "free")), pch = 21,
              ylab = "Class within school", type = c("p", "a"),
              xlab = "Mathematics gain from kindergarten to grade 1",
              jitter.y = TRUE))


###################################################
### chunk number 40: fm5
###################################################
(fm5 <- lmer(mathgain ~ 1 + (1|classid) + (1|schoolid), classroom))


###################################################
### chunk number 41: fm6
###################################################
print(fm6 <- lmer(mathgain ~ 1 + mathkind + minority + sex + ses + housepov
                  + (1|classid) + (1|schoolid), classroom), corr = FALSE)


###################################################
### chunk number 42: fm7
###################################################
print(fm7 <- lmer(mathgain ~ 1 + mathkind + minority + ses + housepov
                  + (1|classid) + (1|schoolid), classroom), corr = FALSE)


###################################################
### chunk number 43: fm8
###################################################
print(fm8 <- lmer(mathgain ~ mathkind + minority + ses
                  + (1|classid) + (1|schoolid), classroom), corr = FALSE)


###################################################
### chunk number 44: Classpredi
###################################################
print(dotplot(ranef(fm8, post = "TRUE"), strip = FALSE,
              scales = list(y = list(draw = FALSE)))$classid)


###################################################
### chunk number 45: classpred2show
###################################################
qqmath(ranef(fm8, post = TRUE))$classid


###################################################
### chunk number 46: Classpred2
###################################################
print(qqmath(ranef(fm8, post = TRUE),strip=FALSE)$classid)


###################################################
### chunk number 47: Schoolpred
###################################################
print(qqmath(ranef(fm8, post = TRUE),strip=FALSE)$schoolid)


###################################################
### chunk number 48: fm9
###################################################
print(fm9M <- lmer(mathgain ~ mathkind + minority + ses
                   + (1|schoolid), classroom, REML = FALSE), corr = FALSE)


###################################################
### chunk number 49: fm8Manova
###################################################
fm8M <- update(fm8, REML = FALSE)
anova(fm9M, fm8M)


