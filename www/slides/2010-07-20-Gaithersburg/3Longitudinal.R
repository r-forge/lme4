###################################################
### chunk number 1: preliminaries
###################################################
options(width=65,show.signif.stars=FALSE)
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
### chunk number 9: fm3anova
###################################################
fm3 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
anova(fm3, fm2)


###################################################
### chunk number 10: rr1
###################################################
(rr2 <- ranef(fm2))


###################################################
### chunk number 11: rr2plot
###################################################
print(plot(rr2, aspect = 1, type = c("g", "p"))[[1]])


###################################################
### chunk number 12: shrinkage
###################################################
df <- coef(lmList(Reaction ~ Days | Subject, sleepstudy))
cc1 <- as.data.frame(coef(fm2)[["Subject"]])
names(cc1) <- c("A", "B")
df <- cbind(df, cc1)
with(df,
     print(xyplot(`(Intercept)` ~ Days, aspect = 1,
            x1 = B, y1 = A, 
            panel = function(x, y, x1, y1, subscripts, ...) {
                panel.grid(h = -1, v = -1)
                x1 <- x1[subscripts]
                y1 <- y1[subscripts]
                panel.arrows(x, y, x1, y1, type = "closed", length = 0.1,
                             angle = 15, ...)
                panel.points(x, y,
                             col = trellis.par.get("superpose.symbol")$col[2])
                panel.points(x1, y1)
            },
                   key = list(space = "top", columns = 2,
                   text = list(c("Mixed model", "Within-group")),
                   points = list(col = trellis.par.get("superpose.symbol")$col[1:2],
                   pch = trellis.par.get("superpose.symbol")$pch[1:2]))
               )))


###################################################
### chunk number 13: shrinkfit
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
             },
             index.cond = function(x,y) coef(lm(y ~ x))[1],
             xlab = "Days of sleep deprivation",
             ylab = "Average reaction time (ms)",
             key = list(space = "top", columns = 2,
             text = list(c("Mixed model", "Within-group")),
             lines = list(col = trellis.par.get("superpose.line")$col[2:1],
             lty = trellis.par.get("superpose.line")$lty[2:1]))))


###################################################
### chunk number 14: caterpillar
###################################################
print(dotplot(ranef(fm1,post=TRUE),
              scales = list(x = list(relation = 'free')))[["Subject"]])


