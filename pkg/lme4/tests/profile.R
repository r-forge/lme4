library(lme4)

fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)

## 0.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
##
system.time( tpr <- profile(fm01ML) )

## test all combinations of 'which'
wlist <- list(1:3,1:2,1,2:3,2,3,c(1,3))
invisible(lapply(wlist,function(w) xyplot(profile(fm01ML,which=w))))
tpr2 <- profile(fm01ML,which=1:3)
tpr3 <- profile(fm01ML,which=2:3)  ## can't plot
tpr4 <- profile(fm01ML,which=3)
tpr5 <- profile(fm01ML,which=2)  ## can't plot
tpr6 <- profile(fm01ML,which=1)

(confint(tpr) -> CIpr)
print(xyplot(tpr))
##  comparing against lme4a reference values -- but lme4 returns sigma
## rather than log(sigma)
stopifnot(dim(CIpr) == c(3,2),
          all.equal(unname(CIpr[".sigma",]),exp(c(3.64362, 4.21446)), tol=1e-6),
          all.equal(unname(CIpr["(Intercept)",]),c(1486.451500,1568.548494)))

## 2D profiles
fm2ML <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin, REML=0)
system.time(pr2 <- profile(fm2ML))
(confint(pr2) -> CIpr2)

lme4a_CIpr2 <-
structure(c(0.633565787613112, 1.09578224011285, -0.721864513060904,
21.2666273835452, 1.1821039843372, 3.55631937954106, -0.462903300019305,
24.6778176174587), .Dim = c(4L, 2L), .Dimnames = list(c(".sig01",
".sig02", ".lsig", "(Intercept)"), c("2.5 %", "97.5 %")))
lme4a_CIpr2[".lsig",] <- exp(lme4a_CIpr2[".lsig",])

stopifnot(all.equal(unname(CIpr2),unname(lme4a_CIpr2),tol=1e-6))

print(xyplot(pr2, absVal=0, aspect=1.3, layout=c(4,1)))
print(splom(pr2))

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
system.time(pr4 <- profile(gm1))  ## ~ 7 seconds

xyplot(pr4,layout=c(5,1),as.table=TRUE)
splom(pr4)

nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
             Orange, start = c(Asym = 200, xmid = 725, scal = 350))
## pr5 <- profile(nm1)
## xyplot(.zeta~.focal|.par,type="l",data=as.data.frame(pr5),
##        scale=list(x=list(relation="free")),
##       as.table=TRUE)


## NOT RUN:  ~ 4 theta-variables, 19 seconds
fm3ML <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML=FALSE)
if (FALSE) {
  system.time(pr3 <- profile(fm3ML))
  xyplot(pr3)
  print(splom(pr3))
}
