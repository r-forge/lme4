###################################################
### chunk number 1: preliminaries
###################################################
#line 17 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
options(width=65,show.signif.stars=FALSE,str=strOptions(strict.width="cut"))
library(Rcpp)
library(minqa)
library(lattice)
library(Matrix)
library(MatrixModels)
library(lme4a)
data(Multilocation, package = "SASmixed")
attr(Multilocation, "ginfo") <-NULL
data(Early, package="mlmRev")
#lattice.options(default.theme = function() standard.theme())
lattice.options(default.theme = function() standard.theme(color=FALSE))
if (file.exists("fm11.rda")) {
    load("fm11.rda")
} else {
    fm11 <- lmer(Adj ~ Trt + (0+Trt|Location) + (1|Grp), Multilocation, REML=FALSE)
    save(fm11, file="fm11.rda")
}


###################################################
### chunk number 2: sleepxy
###################################################
#line 85 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
                    layout = c(9,2), type = c("g", "p", "r"),
                    index.cond = function(x,y) coef(lm(y ~ x))[1],
                    xlab = "Days of sleep deprivation",
                    ylab = "Average reaction time (ms)"))


###################################################
### chunk number 3: Sl1
###################################################
#line 137 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(plot(confint(lmList(Reaction ~ Days | Subject, sleepstudy),
                   pooled = TRUE), order = 1))


###################################################
### chunk number 4: fm1
###################################################
#line 149 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)


###################################################
### chunk number 5: sleepZ
###################################################
#line 167 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(image(fm1@re@Zt,xlab=NULL,ylab=NULL,sub=NULL))


###################################################
### chunk number 6: sm1
###################################################
#line 174 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))


###################################################
### chunk number 7: fm2
###################################################
#line 215 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))


###################################################
### chunk number 8: anovafm1fm2
###################################################
#line 231 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
anova(fm2, fm1)


###################################################
### chunk number 9: rr1
###################################################
#line 276 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
(rr2 <- ranef(fm2))


###################################################
### chunk number 10: rr2plot
###################################################
#line 283 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(plot(rr2, aspect = 1, type = c("g", "p"))[[1]])


###################################################
### chunk number 11: shrinkage
###################################################
#line 310 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
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
#line 352 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
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
#line 381 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(dotplot(ranef(fm1,post=TRUE),
              scales = list(x = list(relation = 'free')))[["Subject"]])


###################################################
### chunk number 14: Multilocation
###################################################
#line 442 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
str(Multilocation)


###################################################
### chunk number 15: Multiplot1
###################################################
#line 450 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(dotplot(reorder(Grp, Adj) ~ Adj, Multilocation,
              groups=Trt, type=c("p","a"),
              auto.key=list(columns=4,lines=TRUE)))


###################################################
### chunk number 16: Multiplot2
###################################################
#line 467 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
ll <- with(Multilocation, reorder(Location, Adj))
print(dotplot(reorder(reorder(Grp, Adj), as.numeric(ll)) ~ Adj|ll, Multilocation,
              groups=Trt, type=c("p","a"), strip=FALSE, strip.left=TRUE, layout=c(1,9),
              auto.key=list(columns=4,lines=TRUE),
              scales = list(y=list(relation="free"))))


###################################################
### chunk number 17: fm3
###################################################
#line 501 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(fm3 <- lmer(Adj ~ Trt + (1|Grp), Multilocation), corr=FALSE)


###################################################
### chunk number 18: lrt
###################################################
#line 519 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
fm4 <- lmer(Adj ~ 1 + (1|Grp), Multilocation)
anova(fm4, fm3)


###################################################
### chunk number 19: fm5
###################################################
#line 536 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
anova(fm5 <- lmer(Adj ~ Location + Trt + (1|Grp), Multilocation))


###################################################
### chunk number 20: fm6
###################################################
#line 548 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
anova(fm6 <- lmer(Adj ~ Location*Trt + (1|Grp), Multilocation))
anova(fm5, fm6)


###################################################
### chunk number 21: fm7
###################################################
#line 557 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(fm7 <- lmer(Adj ~ Trt + (1|Location) + (1|Grp), Multilocation), corr = FALSE)


###################################################
### chunk number 22: fm8
###################################################
#line 570 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
fm8 <- lmer(Adj ~ Trt + (1|Location), Multilocation)
anova(fm8, fm7)


###################################################
### chunk number 23: fm9
###################################################
#line 599 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
(fm9 <- lmer(Adj ~ Trt + (1|Trt:Location) + (1|Location), Multilocation, REML=FALSE))


###################################################
### chunk number 24: fm10
###################################################
#line 607 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
(fm10 <- update(fm9, . ~ . + (1|Grp)))


###################################################
### chunk number 25: anovafm10
###################################################
#line 615 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
anova(fm10, fm8)


###################################################
### chunk number 26: fm11 eval=FALSE
###################################################
## #line 639 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
## fm11 <- lmer(Adj ~ Trt + (Trt|Location) + (1|Grp), Multilocation, REML=FALSE)


###################################################
### chunk number 27: fm11 eval=FALSE
###################################################
## #line 643 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
## fm11 <- lmer(Adj ~ Trt + (0+Trt|Location) + (1|Grp), Multilocation, REML=FALSE)


###################################################
### chunk number 28: 
###################################################
#line 653 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
cat(paste(capture.output(print(fm11))[4:15], collapse="\n"), "\n")


###################################################
### chunk number 29: fm11kappa
###################################################
#line 700 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
kappa(fm11@re@Lambda)
rcond(fm11@re@Lambda)


###################################################
### chunk number 30: fm11verb
###################################################
#line 714 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
fm11 <- lmer(Adj ~ Trt + (0+Trt|Location) + (1|Grp), Multilocation, REML=FALSE, verbose=TRUE)


###################################################
### chunk number 31: fm11
###################################################
#line 717 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
fm11@re@theta


###################################################
### chunk number 32: EarlyData
###################################################
#line 746 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(xyplot(cog ~ age | id, Early, type = c("g",'b'), aspect = 'xy',
             layout = c(29,4), between = list(y = c(0,0.5)),
#             skip = rep(c(FALSE,TRUE),c(58,11)),
             xlab = "Age (yr)",
             ylab = "Cognitive development score",
             scales = list(x = list(tick.number = 3, alternating = TRUE,
                           labels = c("1","","2"), at = c(1,1.5,2))),
             par.strip.text = list(cex = 0.7)))


###################################################
### chunk number 33: fm12
###################################################
#line 766 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
Early <- within(Early, tos <- age-0.5)
fm12 <- lmer(cog ~ tos+trt:tos+(tos|id), Early, verbose=TRUE)


###################################################
### chunk number 34: fm12show
###################################################
#line 775 "/home/bates/Documents/slides/2011-03-16-Amsterdam/2Longitudinal.Rnw"
print(fm12, corr=FALSE)


