###################################################
### chunk number 1: preliminaries
###################################################
options(width=60,show.signif.stars=FALSE)
library(lattice)
library(Matrix)
library(lme4)


###################################################
### chunk number 2: SSmicmenplot
###################################################
  xx <- seq(0, 5, len = 101)
  yy <- 5 * xx/(1+xx)
  par(mar = c(0, 0, 0, 0))
  plot(xx, yy, type = "l", axes = F, ylim = c(0,6), xlim = c(-1, 5),
       xlab = "", ylab = "", lwd = 2)
  usr <- par("usr")
  arrows(usr[1], 0, usr[2], 0, length = 0.1, angle = 25)
  arrows(0, usr[3], 0, usr[4], length = 0.1, angle = 25)
  text(usr[2] - 0.2, 0.1, "x", adj = c(1, 0))
  text(-0.1, usr[4], "y", adj = c(1, 1))
  abline(h = 5, lty = 2, lwd = 0)
  arrows(-0.8, 2.1, -0.8, 0, length = 0.1, angle = 25)
  arrows(-0.8, 2.9, -0.8, 5, length = 0.1, angle = 25)
  text(-0.8, 2.5, expression(phi[1]), adj = c(0.5, 0.5))
  segments(1, 0, 1, 2.7, lty = 2, lwd = 0.75)
  text(1, 2.7, expression(phi[2]), adj = c(0.5, 0))


###################################################
### chunk number 3: ssasympplot
###################################################
xx <- seq(0, 5, len = 101)
yy <- 5 - 4 * exp(-xx/(2*log(2)))
par(mar = c(0, 0, 0, 0))
plot(xx, yy, type = "l", axes = F, ylim = c(0,6), xlim = c(-1, 5),
xlab = "", ylab = "", lwd = 2)
usr <- par("usr")
arrows(usr[1], 0, usr[2], 0, length = 0.1, angle = 25)
arrows(0, usr[3], 0, usr[4], length = 0.1, angle = 25)
text(usr[2] - 0.2, 0.1, "x", adj = c(1, 0))
text(-0.1, usr[4], "y", adj = c(1, 1))
abline(h = 5, lty = 2, lwd = 0)
arrows(-0.8, 2.1, -0.8, 0, length = 0.1, angle = 25)
arrows(-0.8, 2.9, -0.8, 5, length = 0.1, angle = 25)
text(-0.8, 2.5, expression(phi[1]), adj = c(0.5, 0.5))
segments(-0.4, 1, 0, 1, lty = 2, lwd = 0.75)
arrows(-0.3, 0.25, -0.3, 0, length = 0.07, angle = 25)
arrows(-0.3, 0.75, -0.3, 1, length = 0.07, angle = 25)
text(-0.3, 0.5, expression(phi[2]), adj = c(0.5, 0.5))
segments(1, 3.025, 1, 4, lty = 2, lwd = 0.75)
arrows(0.2, 3.5, 0, 3.5, length = 0.08, angle = 25)
arrows(0.8, 3.5, 1, 3.5, length = 0.08, angle = 25)
text(0.5, 3.5, expression(t[0.5]), adj = c(0.5, 0.5))


###################################################
### chunk number 4: SSlogisplot
###################################################
xx <- seq(-0.5, 5, len = 101)
yy <- 5 / ( 1 + exp((2-xx)))
par(mar = c(0, 0, 0, 0))
plot(xx, yy, type = "l", axes = F, ylim = c(0,6), xlim = c(-1, 5),
     xlab = "", ylab = "", lwd = 2)
usr <- par("usr")
arrows(usr[1], 0, usr[2], 0, length = 0.1, angle = 25)
arrows(0, usr[3], 0, usr[4], length = 0.1, angle = 25)
text(usr[2] - 0.2, 0.1, "x", adj = c(1, 0))
text(-0.1, usr[4], "y", adj = c(1, 1))
abline(h = 5, lty = 2, lwd = 0)
arrows(-0.8, 2.1, -0.8, 0, length = 0.1, angle = 25)
arrows(-0.8, 2.9, -0.8, 5, length = 0.1, angle = 25)
text(-0.8, 2.5, expression(phi[1]), adj = c(0.5, 0.5))
segments(2, 0, 2, 4.0, lty = 2, lwd = 0.75)
text(2, 4.0, expression(phi[2]), adj = c(0.5, 0))
segments(3, 5/(1+exp(-1)) + 0.025, 3, 4.0, lty = 2, lwd = 0.75)
arrows(2.3, 3.8, 2.0, 3.8, length = 0.08, angle = 25)
arrows(2.7, 3.8, 3.0, 3.8, length = 0.08, angle = 25)
text(2.5, 3.8, expression(phi[3]), adj = c(0.5, 0.5))


###################################################
### chunk number 5: Orangeplot
###################################################
print(xyplot(circumference ~ age, Orange, groups = Tree, type = c("g","b"),
             xlab = "Age of tree (days)", ylab = "Circumference"))


###################################################
### chunk number 6: nm1
###################################################
print(nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~
                   Asym | Tree, Orange,
                   start = c(Asym = 200, xmid = 770, scal = 120)),
      corr = FALSE)


###################################################
### chunk number 7: OrangeREplot
###################################################
print(dotplot(ranef(nm1, post = TRUE),
              ylab = "Tree", strip = FALSE)[[1]])


###################################################
### chunk number 8: nm2
###################################################
print(nm2 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~
                   (Asym | Tree) + (xmid | Tree) + (scal|Tree), Orange,
                   start = c(Asym = 200, xmid = 770, scal = 120)),
      corr = FALSE)


###################################################
### chunk number 9: nm3
###################################################
print(nm3 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~
                   (Asym + scal|Tree), Orange,
                   start = c(Asym = 200, xmid = 770, scal = 120)),
      corr = FALSE)


###################################################
### chunk number 10: Orangereplot2
###################################################
print(dotplot(ranef(nm3, post = TRUE),
              scales = list(x = list(relation = "free")),
              ylab = "Tree")[[1]])


###################################################
### chunk number 11: Theophplot
###################################################
print(xyplot(conc ~ Time|Subject, Theoph, type = c("g","b"),
             xlab = "Time since drug administration (hr)",
             ylab = "Serum concentration (mg/l)"))


###################################################
### chunk number 12: nm4
###################################################
Th.start <- c(lKe = -2.5, lKa = 0.5, lCl = -3)
print(nm4 <- nlmer(conc ~ SSfol(Dose,Time,lKe,lKa,lCl) ~
                   (lKe+lKa+lCl|Subject), Theoph, start = Th.start),
      corr = FALSE)


###################################################
### chunk number 13: nm5
###################################################
print(nm5 <- nlmer(conc ~ SSfol(Dose,Time,lKe,lKa,lCl) ~
                   (lKa+lCl|Subject), Theoph, start = Th.start),
      corr = FALSE)


###################################################
### chunk number 14: 
###################################################
print(nm6 <- nlmer(conc ~ SSfol(Dose,Time,lKe,lKa,lCl) ~
                   (lKa|Subject) + (lCl|Subject), Theoph,
                   start = Th.start),
      corr = FALSE)


###################################################
### chunk number 15: TheophREplot
###################################################
print(dotplot(ranef(nm5, post = TRUE), ylab = "Subject",
              scales=list(x=list(relation="free")))[[1]])


