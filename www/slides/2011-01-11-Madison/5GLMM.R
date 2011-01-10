###################################################
### chunk number 1: preliminaries
###################################################
#line 17 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
options(width = 70, show.signif.stars = FALSE)
data(Contraception, package = "mlmRev")
library(lattice)
library(Matrix)
library(MatrixModels)
library(Rcpp)
library(minqa)
library(lme4a)
lattice.options(default.theme = function() standard.theme())
#lattice.options(default.theme = function() standard.theme(color=FALSE))


###################################################
### chunk number 2: BernoulliLink
###################################################
#line 218 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
logit <- function(mu) {
  mu <- pmax(.Machine$double.eps, pmin(1-.Machine$double.eps, mu))
  log(mu/(1-mu))
}
mu <- seq(0.001, 0.999, len = 999)
print(xyplot(logit(mu) ~ mu, type = c("g", "l"), 
             xlab = expression(mu), 
             ylab = expression(eta == log(frac(mu, 1-mu)))))


###################################################
### chunk number 3: BernoulliinvLink
###################################################
#line 232 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
linkinv <- function(eta) 1/(1+exp(-eta))
eta <- seq(-7,7,len = 701)
print(xyplot(linkinv(eta) ~ eta, type = c("g","l"), 
             xlab = expression(eta),
             ylab = expression(mu == frac(1,1+exp(-eta)))))


###################################################
### chunk number 4: Contra1
###################################################
#line 445 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(xyplot(ifelse(use == "Y", 1, 0) ~ age|urban, Contraception,
             groups = livch, type = c("g", "smooth"),
             auto.key = list(space = "top", points = FALSE,
             lines = TRUE, columns = 4),
             ylab = "Proportion", xlab = "Centered age"))


###################################################
### chunk number 5: cm1
###################################################
#line 479 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(cm1 <- glmer(use ~ age + I(age^2) + urban + livch + (1|district), 
             Contraception, binomial), corr = FALSE)


###################################################
### chunk number 6: ch
###################################################
#line 503 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
Contraception$ch <- factor(Contraception$livch != 0, labels = c("N","Y"))


###################################################
### chunk number 7: cm2
###################################################
#line 508 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(cm2 <- glmer(use ~ age + I(age^2) + urban + ch + (1|district),
                  Contraception, binomial), corr = FALSE)


###################################################
### chunk number 8: anovac
###################################################
#line 518 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
anova(cm2, cm1)


###################################################
### chunk number 9: Contra2
###################################################
#line 533 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(xyplot(ifelse(use == "Y", 1, 0) ~ age|urban, Contraception,
             groups = ch, type = c("g", "smooth"),
             auto.key = list(space = "top", points = FALSE,
             lines = TRUE, columns = 2),
             ylab = "Proportion", xlab = "Centered age"))


###################################################
### chunk number 10: cm3
###################################################
#line 545 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(cm3 <- glmer(use ~ age*ch + I(age^2) + urban + (1|district),
                   Contraception, binomial), corr = FALSE)


###################################################
### chunk number 11: ContraCat
###################################################
#line 554 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(qqmath(ranef(cm3, post=TRUE), strip=FALSE)[[1]])


###################################################
### chunk number 12: urbanRural
###################################################
#line 566 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
cat(head(capture.output(xtabs(~urban+district, Contraception)),7),sep='\n')


###################################################
### chunk number 13: cm4
###################################################
#line 573 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
(cm4 <- glmer(use ~ age*ch + I(age^2) + urban + (urban|district),
              Contraception, binomial))


###################################################
### chunk number 14: anovacm4
###################################################
#line 581 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
anova(cm4,cm3)


###################################################
### chunk number 15: ContraCat2
###################################################
#line 596 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(qqmath(ranef(cm4, post = TRUE))$district)


###################################################
### chunk number 16: ContraSc
###################################################
#line 605 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(plot(ranef(cm4), type = c("g","p"), aspect = 1)$district)


###################################################
### chunk number 17: cm5
###################################################
#line 612 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(cm5 <- glmer(use ~ age*ch + I(age^2) + urban + (1|urban:district) + (1|district),
                   Contraception, binomial), corr=FALSE)


###################################################
### chunk number 18: cm6
###################################################
#line 620 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
print(cm6 <- glmer(use ~ age*ch + I(age^2) + urban + (1|urban:district),
                   Contraception, binomial), corr=FALSE)


###################################################
### chunk number 19: anovacm654
###################################################
#line 628 "/home/bates/Documents/slides/2011-01-11-Madison/5GLMM.Rnw"
anova(cm6,cm5,cm4)


