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
fm9 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
            sleepstudy, REML = 0)
if (file.exists("pr9.rda")) {
    load("pr9.rda")
} else {
    pr9 <- profile(fm9)
    save(pr9, file = "pr9.rda")
}
data(Orthodont, package = "MEMSS")


###################################################
### chunk number 2: strsleepstudy
###################################################
str(sleepstudy)


###################################################
### chunk number 3: sleepxyplot
###################################################
print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
             layout = c(6,3), type = c("g", "p", "r"),
             index.cond = function(x,y) coef(lm(y ~ x))[1],
             xlab = "Days of sleep deprivation",
             ylab = "Average reaction time (ms)"))


###################################################
### chunk number 4: xtabssleep
###################################################
xtabs(~ Subject + Days, sleepstudy)


###################################################
### chunk number 5: fm8
###################################################
(fm8 <- lmer(Reaction ~ 1 + Days + (1 + Days|Subject), sleepstudy,
             REML = 0))


###################################################
### chunk number 6: modesfm8
###################################################
head(ranef(fm8)[["Subject"]])


###################################################
### chunk number 7: fm8reshow
###################################################
cat(paste(capture.output(print(fm8))[8:11], collapse = "\n"), "\n")


###################################################
### chunk number 8: fm8LambdaL
###################################################
print(image(env(fm8)$Lambda, sub = NULL, xlab = expression(Lambda),
            ylab = NULL),
      split = c(1,1,3,1), more = TRUE)
print(image(tcrossprod(env(fm8)$Lambda), sub = NULL,
            xlab = expression(Sigma), ylab = NULL),
      split = c(2,1,3,1), more = TRUE)
print(image(get("L", envir = env(fm8)),
            sub = NULL, xlab = "L", ylab = NULL),
      split = c(3,1,3,1))


###################################################
### chunk number 9: formulafm9bad
###################################################
Reaction ~ 1 + Days + (1|Subject) + (Days|Subject)


###################################################
### chunk number 10: fm9
###################################################
(fm9 <- lmer(Reaction ~ 1 + Days + (1|Subject) + (0+Days|Subject),
             sleepstudy, REML = 0))


###################################################
### chunk number 11: modesfm9
###################################################
head(ranef(fm9)[["Subject"]])


###################################################
### chunk number 12: fm9reshow
###################################################
cat(paste(capture.output(print(fm9))[8:11], collapse = "\n"), "\n")


###################################################
### chunk number 13: fm9LambdaL
###################################################
print(image(env(fm9)$Lambda, sub = NULL, xlab = expression(Lambda),
            ylab = NULL),
      split = c(1,1,3,1), more = TRUE)
print(image(tcrossprod(env(fm9)$Lambda), sub = NULL,
            xlab = expression(Sigma), ylab = NULL),
      split = c(2,1,3,1), more = TRUE)
print(image(get("L", envir = env(fm9)),
            sub = NULL, xlab = "L", ylab = NULL),
      split = c(3,1,3,1))


###################################################
### chunk number 14: fm89Zt
###################################################
print(image(env(fm8)$Zt, sub = NULL, xlab = NULL, ylab = NULL,
            scales = list(x = list(labels = NULL))),
      pos = c(0,0.38,1,1), more = TRUE)
print(image(env(fm9)$Zt, sub = NULL, xlab = NULL, ylab = NULL,
            scales = list(x = list(labels = NULL))),
      pos = c(0,0,1,0.62))


###################################################
### chunk number 15: formulafm8
###################################################
Reaction ~ Days + (Days|Subject)


###################################################
### chunk number 16: fm1redux
###################################################
Yield ~ (1|Batch)


###################################################
### chunk number 17: fm1redux
###################################################
Yield ~ 1|Batch


###################################################
### chunk number 18: anovafm9fm8
###################################################
anova(fm9, fm8)


###################################################
### chunk number 19: pr9plt
###################################################
print(xyplot(pr9, aspect = 1.3, layout = c(5,1)))


###################################################
### chunk number 20: confintpr9
###################################################
confint(pr9)


###################################################
### chunk number 21: pr9pairs
###################################################
print(splom(pr9))


###################################################
### chunk number 22: rr1
###################################################
str(rr1 <- ranef(fm9))


###################################################
### chunk number 23: ranefcoeffm9
###################################################
print(plot(rr1, type = c("g","p"), aspect = 1)$Subject,
      split = c(1,1,2,1), more = TRUE)
print(plot(coef(fm9), type = c("g","p"), aspect = 1)$Subject,
      split = c(2,1,2,1))


###################################################
### chunk number 24: shrinkage
###################################################
df <- coef(lmList(Reaction ~ Days | Subject, sleepstudy))
fclow <- subset(df, `(Intercept)` < 251)
fchigh <- subset(df, `(Intercept)` > 251)
cc1 <- as.data.frame(coef(fm9)$Subject)
names(cc1) <- c("A", "B")
df <- cbind(df, cc1)
ff <- fixef(fm9)
with(df,
     print(xyplot(`(Intercept)` ~ Days, aspect = 1,
                  x1 = B, y1 = A,
                  panel = function(x, y, x1, y1, subscripts, ...) {
                      panel.grid(h = -1, v = -1)
                      x1 <- x1[subscripts]
                      y1 <- y1[subscripts]
                      larrows(x, y, x1, y1, type = "closed", length = 0.1,
                              angle = 15, ...)
                      lpoints(x, y,
                              pch = trellis.par.get("superpose.symbol")$pch[2],
                              col = trellis.par.get("superpose.symbol")$col[2])
                      lpoints(x1, y1,
                              pch = trellis.par.get("superpose.symbol")$pch[1],
                              col = trellis.par.get("superpose.symbol")$col[1])
                      lpoints(ff[2], ff[1], 
                              pch = trellis.par.get("superpose.symbol")$pch[3],
                              col = trellis.par.get("superpose.symbol")$col[3])
                      ltext(fclow[,2], fclow[,1], row.names(fclow),
                            adj = c(0.5, 1.7))
                      ltext(fchigh[,2], fchigh[,1], row.names(fchigh),
                            adj = c(0.5, -0.6))
                  },
                  key = list(space = "top", columns = 3,
                  text = list(c("Mixed model", "Within-group", "Population")),
                   points = list(col = trellis.par.get("superpose.symbol")$col[1:3],
                   pch = trellis.par.get("superpose.symbol")$pch[1:3]))
               )))


###################################################
### chunk number 25: shrinkfit
###################################################
print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
             layout = c(6,3), type = c("g", "p", "r"),
             coef.list = df[,3:4],
             panel = function(..., coef.list) {
                 panel.xyplot(...)
                 panel.abline(as.numeric(coef.list[packet.number(),]),
                              col.line = trellis.par.get("superpose.line")$col[2],
                              lty = trellis.par.get("superpose.line")$lty[2]
                              )
                 panel.abline(fixef(fm9),
                              col.line = trellis.par.get("superpose.line")$col[4],
                              lty = trellis.par.get("superpose.line")$lty[4]
                              )
             },
             index.cond = function(x,y) coef(lm(y ~ x))[1],
             xlab = "Days of sleep deprivation",
             ylab = "Average reaction time (ms)",
             key = list(space = "top", columns = 3,
             text = list(c("Within-subject", "Mixed model", "Population")),
             lines = list(col = trellis.par.get("superpose.line")$col[c(1:2,4)],
             lty = trellis.par.get("superpose.line")$lty[c(1:2,4)]))))


###################################################
### chunk number 26: caterpillar
###################################################
print(dotplot(ranef(fm8,post=TRUE),
              scales = list(x = list(relation = "free")))[[1]])


###################################################
### chunk number 27: orthofem
###################################################
print(xyplot(distance ~ age|Subject, Orthodont, subset = Sex == "Female",
             index.cond = function(x,y) y[x == 8],
             aspect = 'xy', layout = c(11,1), type = c("g","p","r"),
             xlab = "Age (yr)",
             ylab = "Distance from pituitary to pterygomaxillary fissure (mm)"))


