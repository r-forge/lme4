###################################################
### chunk number 1: preliminaries
###################################################
options(width=65,show.signif.stars=FALSE,str=strOptions(strict.width="cut"))
library(lattice)
library(Matrix)
library(lme4a)
lattice.options(default.theme = function() standard.theme())
#lattice.options(default.theme = function() standard.theme(color=FALSE))


###################################################
### chunk number 2: sleepxy
###################################################
print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
                    layout = c(9,2), type = c("g", "p", "r"),
                    index.cond = function(x,y) coef(lm(y ~ x))[1],
                    xlab = "Days of sleep deprivation",
                    ylab = "Average reaction time (ms)"))


###################################################
### chunk number 3: Sl1
###################################################
print(plot(confint(lmList(Reaction ~ Days | Subject, sleepstudy),
                   pooled = TRUE), order = 1))


###################################################
### chunk number 4: fm1
###################################################
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)


###################################################
### chunk number 5: sleepZ
###################################################
print(image(fm1@re@Zt,xlab=NULL,ylab=NULL,sub=NULL))


###################################################
### chunk number 6: sm1
###################################################
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))


###################################################
### chunk number 7: fm2
###################################################
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))


###################################################
### chunk number 8: anovafm1fm2
###################################################
anova(fm2, fm1)


###################################################
### chunk number 9: rr1
###################################################
(rr2 <- ranef(fm2))


###################################################
### chunk number 10: rr2plot
###################################################
print(plot(rr2, aspect = 1, type = c("g", "p"))[[1]])


###################################################
### chunk number 11: shrinkage
###################################################
df <- coef(lmList(Reaction ~ Days | Subject, sleepstudy))
fclow <- subset(df, `(Intercept)` < 251)
fchigh <- subset(df, `(Intercept)` > 251)
cc1 <- as.data.frame(coef(fm2)$Subject)
names(cc1) <- c("A", "B")
df <- cbind(df, cc1)
ff <- fixef(fm2)
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
### chunk number 12: shrinkfit
###################################################
print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
             layout = c(9,2), type = c("g", "p", "r"),
             coef.list = df[,3:4],
             panel = function(..., coef.list) {
                 panel.xyplot(...)
                 panel.abline(as.numeric(coef.list[packet.number(),]),
                              col.line = trellis.par.get("superpose.line")$col[2],
                              lty = trellis.par.get("superpose.line")$lty[2]
                              )
                 panel.abline(fixef(fm2),
                              col.line = trellis.par.get("superpose.line")$col[4],
                              lty = trellis.par.get("superpose.line")$lty[4]
                              )
             },
             index.cond = function(x,y) coef(lm(y ~ x))[1],
             xlab = "Days of sleep deprivation",
             ylab = "Average reaction time (ms)",
             key = list(space = "top", columns = 3,
             text = list(c("Within-subject", "Mixed model", "Population")),
             lines = list(col = trellis.par.get("superpose.line")$col[c(2:1,4)],
             lty = trellis.par.get("superpose.line")$lty[c(2:1,4)]))))


###################################################
### chunk number 13: caterpillar
###################################################
print(dotplot(ranef(fm1,post=TRUE),
              scales = list(x = list(relation = 'free')))[["Subject"]])


