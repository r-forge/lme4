###################################################
### chunk number 1: preliminaries
###################################################
options(width = 70, show.signif.stars = FALSE)
data(Contraception, package = "mlmRev")
library(lattice)
library(Matrix)
library(MatrixModels)
library(Rcpp)
library(minqa)
library(lme4a)
lattice.options(default.theme = function() standard.theme())


###################################################
### chunk number 2: Theophplot
###################################################
print(xyplot(conc ~ Time|Subject, Theoph, type = c("g","b"),
             xlab = "Time since drug administration (hr)",
             ylab = "Serum concentration (mg/l)"))


###################################################
### chunk number 3: nm1
###################################################
Th.start <- c(lKe = -2.5, lKa = 0.5, lCl = -3)
nm1 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
             0+lKe+lKa+lCl+(0+lKe|Subject)+(0+lKa|Subject)
             +(0+lCl|Subject), nAGQ=0, Theoph,
             start = Th.start, verbose=TRUE)


###################################################
### chunk number 4: nm1out
###################################################
print(nm1, corr=FALSE)


###################################################
### chunk number 5: nm2
###################################################
nm2 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
             0+lKe+lKa+lCl+(0+lKa|Subject)
             +(0+lCl|Subject), Theoph, nAGQ=0,
             start = Th.start, verbose=TRUE)


###################################################
### chunk number 6: nm3
###################################################
nm3 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
             0+lKe+lKa+lCl+(0+lKa+lCl|Subject),
             Theoph, start = Th.start, verbose=TRUE)


###################################################
### chunk number 7: nm3show
###################################################
print(nm3, corr=FALSE)


###################################################
### chunk number 8: anovanm2nm3
###################################################
anova(nm2,nm3)


