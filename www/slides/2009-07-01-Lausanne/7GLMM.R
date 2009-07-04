>###################################################
### chunk number 1: preliminaries
###################################################
options(width = 70, show.signif.stars = FALSE)
data(Contraception, package = "mlmRev")
library(lattice)
library(Matrix)
library(lme4)
lattice.options(default.theme = function() standard.theme())


###################################################
### chunk number 2: BernoulliLink
###################################################
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
linkinv <- function(eta) 1/(1+exp(-eta))
eta <- seq(-7,7,len = 701)
print(xyplot(linkinv(eta) ~ eta, type = c("g","l"), 
             xlab = expression(eta),
             ylab = expression(mu == frac(1,1+exp(-eta)))))


###################################################
### chunk number 4: Contra1
###################################################
print(xyplot(ifelse(use == "Y", 1, 0) ~ age|urban, Contraception,
             groups = livch, type = c("g", "smooth"),
             auto.key = list(space = "top", points = FALSE,
             lines = TRUE, columns = 4),
             ylab = "Proportion", xlab = "Centered age"))


###################################################
### chunk number 5: cm1
###################################################
print(cm1 <- lmer(use ~ age + I(age^2) + urban + livch + (1|district), 
             Contraception, binomial), corr = FALSE)


###################################################
### chunk number 6: ch
###################################################
Contraception$ch <- factor(Contraception$livch != 0, labels = c("N","Y"))


###################################################
### chunk number 7: cm2
###################################################
print(cm2 <- lmer(use ~ age + I(age^2) + urban + ch + (1|district),
                  Contraception, binomial), corr = FALSE)


###################################################
### chunk number 8: anovac
###################################################
anova(cm2, cm1)


###################################################
### chunk number 9: Contra2
###################################################
print(xyplot(ifelse(use == "Y", 1, 0) ~ age|urban, Contraception,
             groups = ch, type = c("g", "smooth"),
             auto.key = list(space = "top", points = FALSE,
             lines = TRUE, columns = 2),
             ylab = "Proportion", xlab = "Centered age"))


###################################################
### chunk number 10: cm3
###################################################
print(cm3 <- glmer(use ~ age*ch + I(age^2) + urban + (1|district),
                   Contraception, binomial), corr = FALSE)


###################################################
### chunk number 11: ContraCat
###################################################
print(qqmath(ranef(cm3, post = TRUE))[[1]])


###################################################
### chunk number 12: urbanRural
###################################################
cat(head(capture.output(xtabs(~urban+district, Contraception)),7),sep='\n')


###################################################
### chunk number 13: cm4
###################################################
(cm4 <- glmer(use ~ age*ch + I(age^2) + urban + (urban|district),
              Contraception, binomial))


###################################################
### chunk number 14: anovacm4
###################################################
anova(cm4,cm3)


###################################################
### chunk number 15: ContraCat2
###################################################
print(qqmath(ranef(cm4, post = TRUE))$district)


###################################################
### chunk number 16: ContraSc
###################################################
print(plot(ranef(cm4), type = c("g","p"), aspect = 1)$district)


