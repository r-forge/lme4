###################################################
### chunk number 1: preliminaries
###################################################
#line 17 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
options(width=69,show.signif.stars=FALSE,str=strOptions(strict.width="cut"))
library(lattice)
library(Rcpp)
library(minqa)
library(Matrix)
library(MatrixModels)
library(lme4a)
lattice.options(default.theme = function() standard.theme())
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
## #line 57 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
## install.packages("lme4a", repos="http://r-forge.r-project.org")


###################################################
### chunk number 3: require eval=FALSE
###################################################
## #line 63 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
## require(lme4a)


###################################################
### chunk number 4: attach eval=FALSE
###################################################
## #line 67 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
## library(lme4a)


###################################################
### chunk number 5: Dyestuffstr
###################################################
#line 137 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
str(Dyestuff)
summary(Dyestuff)


###################################################
### chunk number 6: Dyestuffplot
###################################################
#line 165 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
set.seed(1234543)
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff,
              ylab = "Batch", jitter.y = TRUE, pch = 21, aspect = 0.32,
              xlab = "Yield of dyestuff (grams of standard color)",
              type = c("p", "a")))


###################################################
### chunk number 7: fm1
###################################################
#line 184 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
fm1 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff)
print(fm1)


###################################################
### chunk number 8: op
###################################################
#line 195 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
op <- options(digits=5)


###################################################
### chunk number 9: extractors
###################################################
#line 206 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
fixef(fm1)
ranef(fm1, drop = TRUE)
fitted(fm1)


###################################################
### chunk number 10: unop
###################################################
#line 213 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
options(op)


###################################################
### chunk number 11: strfm1
###################################################
#line 313 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
str(as.matrix(model.matrix(fm1)))
fm1@re@Zt


###################################################
### chunk number 12: fm1ranef
###################################################
#line 353 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(dotplot(ranef(fm1, postVar = TRUE), strip = FALSE)[[1]])


###################################################
### chunk number 13: update
###################################################
#line 384 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
(fm1M <- update(fm1, REML = FALSE))


###################################################
### chunk number 14: fm1refit
###################################################
#line 413 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
invisible(update(fm1, verbose = TRUE))


###################################################
### chunk number 15: Dyestuff2
###################################################
#line 437 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
str(Dyestuff2)


###################################################
### chunk number 16: Dyestuff2plot
###################################################
#line 444 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff2,
              ylab = "Batch", jitter.y = TRUE, pch = 21, aspect = 0.42,
              xlab = "Simulated response (dimensionless)",
              type = c("p", "a")))


###################################################
### chunk number 17: fm1A
###################################################
#line 458 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
(fm1A <- lmer(Yield ~ 1 + (1|Batch), Dyestuff2, REML=FALSE))


###################################################
### chunk number 18: lm1
###################################################
#line 469 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
summary(lm1 <- lm(Yield ~ 1, Dyestuff2))
logLik(lm1)


###################################################
### chunk number 19: fm1call
###################################################
#line 478 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
fm1@call


###################################################
### chunk number 20: Penicillinstr
###################################################
#line 497 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
str(Penicillin)
xtabs(~ sample + plate, Penicillin)


###################################################
### chunk number 21: PenicillinPlot
###################################################
#line 512 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(dotplot(reorder(plate, diameter) ~ diameter, Penicillin, groups = sample,
              ylab = "Plate", xlab = "Diameter of growth inhibition zone (mm)",
              type = c("p", "a"), auto.key = list(columns = 6, lines = TRUE)))


###################################################
### chunk number 22: fm2
###################################################
#line 522 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
(fm2 <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin))


###################################################
### chunk number 23: 
###################################################
#line 527 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
op <- options(digits = 5)


###################################################
### chunk number 24: fixef2
###################################################
#line 539 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
fixef(fm2)


###################################################
### chunk number 25: ranef2
###################################################
#line 543 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
ranef(fm2, drop = TRUE)


###################################################
### chunk number 26: 
###################################################
#line 548 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
options(op)


###################################################
### chunk number 27: fm2ranef
###################################################
#line 554 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
qrr2 <- dotplot(ranef(fm2, postVar = TRUE), strip = FALSE)
print(qrr2[[1]], pos = c(0,0,1,0.75), more = TRUE)
print(qrr2[[2]], pos = c(0,0.65,1,1))


###################################################
### chunk number 28: fm2Z
###################################################
#line 572 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(image(fm2@re@Zt, xlab = NULL, ylab = NULL, sub = "Z'"))


###################################################
### chunk number 29: fm2formula
###################################################
#line 625 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
fm2@call[["formula"]]


###################################################
### chunk number 30: Pastesstr
###################################################
#line 647 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
str(Pastes)
xtabs(~ batch + sample, Pastes, sparse = TRUE)


###################################################
### chunk number 31: samplegen eval=FALSE
###################################################
## #line 672 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
## Pastes <- within(Pastes, sample <- factor(batch:cask))


###################################################
### chunk number 32: Pastesplot
###################################################
#line 700 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
Pastes <- within(Pastes, bb <- reorder(batch, strength))
Pastes <- within(Pastes, ss <- reorder(reorder(sample, strength),
          as.numeric(batch)))
print(dotplot(ss ~ strength | bb, Pastes,
              strip = FALSE, strip.left = TRUE, layout = c(1, 10),
              scales = list(y = list(relation = "free")),
              ylab = "Sample within batch", type = c("p", "a"),
              xlab = "Paste strength", jitter.y = TRUE))


###################################################
### chunk number 33: fm3
###################################################
#line 714 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
(fm3 <- lmer(strength ~ 1 + (1|batch) + (1|sample), Pastes))


###################################################
### chunk number 34: fm3ranef
###################################################
#line 722 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
qrr3 <- dotplot(ranef(fm3, postVar = TRUE), strip = FALSE)
print(qrr3[[1]], pos = c(0,0,1,0.75), more = TRUE)
print(qrr3[[2]], pos = c(0,0.65,1,1))


###################################################
### chunk number 35: Z3fig
###################################################
#line 743 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(image(fm3@re@Zt, xlab = NULL, ylab = NULL, sub = NULL))


###################################################
### chunk number 36: fm3LRT
###################################################
#line 768 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
fm3M <- update(fm3, REML = FALSE)
fm4M <- lmer(strength ~ 1 + (1|sample),
             Pastes, REML = FALSE)
anova(fm4M, fm3M)


###################################################
### chunk number 37: fm4
###################################################
#line 807 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
(fm4 <- update(fm4M, REML = TRUE))


###################################################
### chunk number 38: xtabsclass
###################################################
#line 870 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
xtabs( ~ xtabs(~ classid, classroom))


###################################################
### chunk number 39: xtabsclass
###################################################
#line 875 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
table(xtabs(~ schoolid,
            unique(subset(classroom, select = c(classid, schoolid)))))


###################################################
### chunk number 40: Schoolsplot
###################################################
#line 884 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
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
### chunk number 41: fm5
###################################################
#line 905 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
(fm5 <- lmer(mathgain ~ 1 + (1|classid) + (1|schoolid), classroom))


###################################################
### chunk number 42: fm6
###################################################
#line 951 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(fm6 <- lmer(mathgain ~ 1 + mathkind + minority + sex + ses + housepov
                  + (1|classid) + (1|schoolid), classroom), corr = FALSE)


###################################################
### chunk number 43: fm7
###################################################
#line 983 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(fm7 <- lmer(mathgain ~ 1 + mathkind + minority + ses + housepov
                  + (1|classid) + (1|schoolid), classroom), corr = FALSE)


###################################################
### chunk number 44: fm8
###################################################
#line 991 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(fm8 <- lmer(mathgain ~ mathkind + minority + ses
                  + (1|classid) + (1|schoolid), classroom), corr = FALSE)


###################################################
### chunk number 45: Classpredi
###################################################
#line 1000 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(dotplot(ranef(fm8, post = "TRUE"), strip = FALSE,
              scales = list(y = list(draw = FALSE)))$classid)


###################################################
### chunk number 46: classpred2show
###################################################
#line 1011 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
qqmath(ranef(fm8, post = TRUE))$classid


###################################################
### chunk number 47: Classpred2
###################################################
#line 1015 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(qqmath(ranef(fm8, post = TRUE),strip=FALSE)$classid)


###################################################
### chunk number 48: Schoolpred
###################################################
#line 1025 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(qqmath(ranef(fm8, post = TRUE),strip=FALSE)$schoolid)


###################################################
### chunk number 49: fm9
###################################################
#line 1035 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
print(fm9M <- lmer(mathgain ~ mathkind + minority + ses
                   + (1|schoolid), classroom, REML = FALSE), corr = FALSE)


###################################################
### chunk number 50: fm8Manova
###################################################
#line 1043 "/home/bates/Documents/slides/2011-01-11-Madison/1Simple.Rnw"
fm8M <- update(fm8, REML = FALSE)
anova(fm9M, fm8M)


