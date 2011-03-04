###################################################
### chunk number 1: preliminaries
###################################################
#line 17 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
options(width=69,show.signif.stars=FALSE,str=strOptions(strict.width="cut"))
library(lattice)
library(Rcpp)
library(minqa)
library(Matrix)
library(MatrixModels)
library(lme4a)
lattice.options(default.theme = function() standard.theme())
#lattice.options(default.theme = function() standard.theme(color=FALSE))
if (file.exists("classroom.rda")) {
    load("classroom.rda")
} else {
    classroom <- within(read.csv("http://www-personal.umich.edu/~bwest/classroom.csv"),
                    {
                        classid <- factor(classid)
                        schoolid <- factor(schoolid)
                        sex <- factor(sex, labels = c("M","F"))
                        minority <- factor(minority, labels = c("N", "Y"))
                    })
    save(classroom, file="classroom.rda")
}


###################################################
### chunk number 2: installlme4 eval=FALSE
###################################################
## #line 58 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
## install.packages("lme4a", repos="http://r-forge.r-project.org/")


###################################################
### chunk number 3: require eval=FALSE
###################################################
## #line 64 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
## require(lme4a)


###################################################
### chunk number 4: attach eval=FALSE
###################################################
## #line 68 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
## library(lme4a)


###################################################
### chunk number 5: Dyestuffstr
###################################################
#line 138 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
str(Dyestuff)
summary(Dyestuff)


###################################################
### chunk number 6: Dyestuffplot
###################################################
#line 166 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
set.seed(1234543)
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff,
              ylab = "Batch", jitter.y = TRUE, pch = 21, aspect = 0.32,
              xlab = "Yield of dyestuff (grams of standard color)",
              type = c("p", "a")))


###################################################
### chunk number 7: fm1
###################################################
#line 185 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
fm1 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff)
print(fm1)


###################################################
### chunk number 8: op
###################################################
#line 196 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
op <- options(digits=5)


###################################################
### chunk number 9: extractors
###################################################
#line 207 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
fixef(fm1)
ranef(fm1, drop = TRUE)
fitted(fm1)


###################################################
### chunk number 10: unop
###################################################
#line 214 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
options(op)


###################################################
### chunk number 11: strfm1
###################################################
#line 314 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
str(as.matrix(model.matrix(fm1)))
fm1@re@Zt


###################################################
### chunk number 12: fm1ranef
###################################################
#line 354 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(dotplot(ranef(fm1, postVar = TRUE), strip = FALSE)[[1]])


###################################################
### chunk number 13: update
###################################################
#line 385 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
(fm1M <- update(fm1, REML = FALSE))


###################################################
### chunk number 14: fm1refit
###################################################
#line 414 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
invisible(update(fm1, verbose = TRUE))


###################################################
### chunk number 15: Dyestuff2
###################################################
#line 438 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
str(Dyestuff2)


###################################################
### chunk number 16: Dyestuff2plot
###################################################
#line 445 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff2,
              ylab = "Batch", jitter.y = TRUE, pch = 21, aspect = 0.42,
              xlab = "Simulated response (dimensionless)",
              type = c("p", "a")))


###################################################
### chunk number 17: fm1A
###################################################
#line 459 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
(fm1A <- lmer(Yield ~ 1 + (1|Batch), Dyestuff2, REML=FALSE))


###################################################
### chunk number 18: lm1
###################################################
#line 470 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
summary(lm1 <- lm(Yield ~ 1, Dyestuff2))
logLik(lm1)


###################################################
### chunk number 19: fm1call
###################################################
#line 479 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
fm1@call


###################################################
### chunk number 20: Penicillinstr
###################################################
#line 498 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
str(Penicillin)
xtabs(~ sample + plate, Penicillin)


###################################################
### chunk number 21: PenicillinPlot
###################################################
#line 513 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(dotplot(reorder(plate, diameter) ~ diameter, Penicillin, groups = sample,
              ylab = "Plate", xlab = "Diameter of growth inhibition zone (mm)",
              type = c("p", "a"), auto.key = list(columns = 6, lines = TRUE)))


###################################################
### chunk number 22: fm2
###################################################
#line 523 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
(fm2 <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin))


###################################################
### chunk number 23: 
###################################################
#line 528 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
op <- options(digits = 5)


###################################################
### chunk number 24: ranef2
###################################################
#line 539 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
ranef(fm2, drop = TRUE)


###################################################
### chunk number 25: 
###################################################
#line 544 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
options(op)


###################################################
### chunk number 26: fm2ranef
###################################################
#line 550 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
qrr2 <- dotplot(ranef(fm2, postVar = TRUE), strip = FALSE)
print(qrr2[[1]], pos = c(0,0,1,0.75), more = TRUE)
print(qrr2[[2]], pos = c(0,0.65,1,1))


###################################################
### chunk number 27: fm2Z
###################################################
#line 568 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(image(fm2@re@Zt, xlab = NULL, ylab = NULL, sub = "Z'"))


###################################################
### chunk number 28: fm2formula
###################################################
#line 621 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
fm2@call[["formula"]]


###################################################
### chunk number 29: Pastesstr
###################################################
#line 643 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
str(Pastes)
xtabs(~ batch + sample, Pastes, sparse = TRUE)


###################################################
### chunk number 30: samplegen eval=FALSE
###################################################
## #line 668 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
## Pastes <- within(Pastes, sample <- factor(batch:cask))


###################################################
### chunk number 31: Pastesplot
###################################################
#line 696 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
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
#line 710 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
(fm3 <- lmer(strength ~ 1 + (1|batch) + (1|sample), Pastes))


###################################################
### chunk number 33: fm3ranef
###################################################
#line 718 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
qrr3 <- dotplot(ranef(fm3, postVar = TRUE), strip = FALSE)
print(qrr3[[1]], pos = c(0,0,1,0.75), more = TRUE)
print(qrr3[[2]], pos = c(0,0.65,1,1))


###################################################
### chunk number 34: Z3fig
###################################################
#line 739 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(image(fm3@re@Zt, xlab = NULL, ylab = NULL, sub = NULL))


###################################################
### chunk number 35: fm3LRT
###################################################
#line 764 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
fm3M <- update(fm3, REML = FALSE)
fm4M <- lmer(strength ~ 1 + (1|sample),
             Pastes, REML = FALSE)
anova(fm4M, fm3M)


###################################################
### chunk number 36: fm4
###################################################
#line 803 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
(fm4 <- update(fm4M, REML = TRUE))


###################################################
### chunk number 37: xtabsclass
###################################################
#line 866 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
xtabs( ~ xtabs(~ classid, classroom))


###################################################
### chunk number 38: xtabsclass
###################################################
#line 871 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
table(xtabs(~ schoolid,
    unique(subset(classroom, select = c(classid, schoolid)))))


###################################################
### chunk number 39: Schoolsplot
###################################################
#line 880 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
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
#line 901 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
(fm5 <- lmer(mathgain ~ 1 + (1|classid) + (1|schoolid),
             classroom))


###################################################
### chunk number 41: fm6
###################################################
#line 948 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(fm6 <- lmer(mathgain ~ 1 + mathkind + minority + sex + ses + housepov
                  + (1|classid) + (1|schoolid), classroom), corr = FALSE)


###################################################
### chunk number 42: fm7
###################################################
#line 980 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(fm7 <- lmer(mathgain ~ 1 + mathkind + minority + ses + housepov
                  + (1|classid) + (1|schoolid), classroom), corr = FALSE)


###################################################
### chunk number 43: fm8
###################################################
#line 988 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(fm8 <- lmer(mathgain ~ mathkind + minority + ses
                  + (1|classid) + (1|schoolid), classroom), corr = FALSE)


###################################################
### chunk number 44: Classpredi
###################################################
#line 997 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(dotplot(ranef(fm8, post = "TRUE"), strip = FALSE,
              scales = list(y = list(draw = FALSE)))$classid)


###################################################
### chunk number 45: classpred2show
###################################################
#line 1008 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
qqmath(ranef(fm8, post = TRUE))$classid


###################################################
### chunk number 46: Classpred2
###################################################
#line 1012 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(qqmath(ranef(fm8, post = TRUE),strip=FALSE)$classid)


###################################################
### chunk number 47: Schoolpred
###################################################
#line 1022 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(qqmath(ranef(fm8, post = TRUE),strip=FALSE)$schoolid)


###################################################
### chunk number 48: fm9
###################################################
#line 1032 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
print(fm9M <- lmer(mathgain ~ mathkind + minority + ses
                   + (1|schoolid), classroom, REML = FALSE), corr = FALSE)


###################################################
### chunk number 49: fm8Manova
###################################################
#line 1040 "/home/bates/Documents/slides/2011-03-16-Amsterdam/1Simple.Rnw"
fm8M <- update(fm8, REML = FALSE)
anova(fm9M, fm8M)


