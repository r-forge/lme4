library(lme4a)
options(show.signif.stars = FALSE)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm1a <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE))
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))

## transformed vars should work[even if non-sensical as here;failed in 0.995-1]
fm2l <- lmer(log(Reaction) ~ log(Days+1) + (log(Days+1)|Subject),
             data = sleepstudy, REML = FALSE)
if(FALSE)# no need for an expand method now
xfm2 <- expand(fm2l)
stopifnot(is(fm1, "merenv"), is(fm2l, "merenv"),
          dim(ranef(fm2l)[[1]]) == c(18, 2),
          TRUE)

## generalized linear mixed model
## (m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##             family = binomial, data = cbpp))
## warnings() ## << FIXME
if(FALSE)#not yet
stopifnot(is(m1,"merenv"), is((cm1 <- coef(m1)), "coef.mer"),
	  dim(cm1$herd) == c(15,4),
          TRUE ## FIXME -- not at all :
	  ## all.equal(fixef(m1),
	  ##           c(-1.39853504914, -0.992334711,
	  ##             -1.12867541477, -1.58037390498), check.attr=FALSE)
	  )

## Simple example by Andrew Gelman (2006-01-10) ----
n.groups <- 10 ; n.reps <- 2
n <- length(group.id <- gl(n.groups, n.reps))
## simulate the varying parameters and the data:
set.seed(0)
a.group <- rnorm(n.groups, 1, 2)
y <- rnorm (n, a.group[group.id], 1)
## fit and summarize the model
fit.1 <- lmer (y ~ 1 + (1 | group.id))
coef (fit.1)# failed in Matrix 0.99-6
(sf1 <- summary(fit.1)) # show() is as without summary()
stopifnot(all.equal(fixef(fit.1), c("(Intercept)" = 1.571312129)),
	  all.equal(ranef(fit.1)[["group.id"]][,"(Intercept)"],
		   c(1.80469, -1.80977, 1.61465, 1.54083, -0.1332,
		     -3.33067, -1.82593, -0.873515, -0.359131, 3.37204),
		    tol = 1e-4)
	  )


## ranef and coef
rr <- ranef(fm1)
stopifnot(is.list(rr), length(rr) == 1, class(rr[[1]]) == "data.frame")
print(plot(rr))
cc <- coef(fm1)
stopifnot(is.list(cc), length(cc) == 1, class(cc[[1]]) == "data.frame")
print(plot(cc))
rr <- ranef(fm2)
stopifnot(is.list(rr), length(rr) == 1, class(rr[[1]]) == "data.frame")
print(plot(rr))
cc <- coef(fm2)
stopifnot(is.list(cc), length(cc) == 1, class(cc[[1]]) == "data.frame")
print(plot(cc))

if (require('MASS', quietly = TRUE)) {
    bacteria$wk2 <- bacteria$week > 2
    contrasts(bacteria$trt) <-
        structure(contr.sdif(3),
                  dimnames = list(NULL, c("diag", "encourage")))
##    print(fm5 <- glmer(y ~ trt + wk2 + (1|ID), bacteria, binomial))
###? the same?
###? print(fm6 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial))

#   if(FALSE) ## FIXME -- binomial  ? -- why the difference?
       ## numbers from 'lme4' ("old"):
#    stopifnot(all.equal(logLik(fm5),
#                        structure(c(ML = -96.13069), nobs = c(n = 220), nall = c(n = 220),
#                                  df = c(p = 5), REML = FALSE, class = "logLik")),
#              all.equal(fixef(fm5),
#			c("(Intercept)"= 2.831609490, "trtdiag"= -1.366722631,
#			  "trtencourage"=0.5840147802, "wk2TRUE"=-1.598591346)))
}

## Invalid factor specification -- used to seg.fault:
set.seed(1)
dat <- data.frame(y = round(10*rnorm(100)), lagoon = factor(rep(1:4,each = 25)),
                  habitat = factor(rep(1:20, each = 5)))
r1  <- lmer(y ~ habitat + (1|habitat:lagoon), data = dat) # ok
stopifnot(all.equal(unname(fixef(r1)),
		    c(1.4, 0, -1, 3.2, -0.6, -5.2, -1, 0.6, -1.2, 1.2,
		      -0.6, 0.8, 3.4, 3.2, -5, -2.6, -2.2, 2.2, 7, -7.4)))

if (FALSE) {   # back to segfaulting again  ----- FIXME !!!!
    try(
        reg <- lmer(y ~ habitat + (1|habitat*lagoon), data = dat) # did seg.fault
        ) # now gives error                 ^- should be ":"
}

## Failure to specify a random effects term - used to give an obscure message
##try(
##m2 <- glmer(incidence / size ~ period, weights = size,
##            family = binomial, data = cbpp)
##)

### mcmcsamp() :
## From: Andrew Gelman <gelman@stat.columbia.edu>
## Date: Wed, 18 Jan 2006 22:00:53 -0500

if (FALSE) {  # mcmcsamp still needs work
    has.coda <- require(coda)
    if(!has.coda)
        cat("'coda' package not available; some outputs will look suboptimal\n")

    ## Very simple example
    y <- 1:10
    group <- gl(2,5)
    (M1 <- lmer (y ~ 1 + (1 | group))) # works fine
    (r1 <- mcmcsamp (M1))              # dito
    r2 <- mcmcsamp (M1, saveb = TRUE)  # gave error in 0.99-* and 0.995-[12]
    (r10 <- mcmcsamp (M1, n = 10, saveb = TRUE))

    ## another one, still simple
    y <- (1:20)*pi
    x <- (1:20)^2
    group <- gl(2,10)
    M1 <- lmer (y ~ 1 | group)
    mcmcsamp (M1, n = 2, saveb=TRUE) # fine

    M2 <- lmer (y ~ 1 + x + (1 + x | group)) # false convergence
    ## should be identical (and is)
    M2 <- lmer (y ~ x + ( x | group))#  false convergence -> simulation doesn't work:
    if(FALSE) ## try(..) fails here (in R CMD check) [[why ??]]
        mcmcsamp (M2, saveb=TRUE)
    ## Error: inconsistent degrees of freedom and dimension ...

    ## mcmc for glmer:
    rG1k <- mcmcsamp(m1, n = 1000)
    summary(rG1k)
    rG2 <- mcmcsamp(m1, n = 3, verbose = TRUE)
}

## Spencer Graves' example (from a post to S-news, 2006-08-03): ----------------
## FIXME?
tstDF <- data.frame(group = letters[1:5], y = 1:5)
var(tstDF$y) # == 2.5
## Now throws an error
try(f.oops <- lmer(y ~ 1 + (1|group), data = tstDF))
##  summary(f.oops) ## or print(Matrix:::formatVC(VarCorr(f.oops)), quote = FALSE)
## ...
##   Groups   Name        Variance Std.Dev.
##   group    (Intercept) 1.81818  1.34840
##   Residual             0.68182  0.82572
## ...
##
##SG>	 This is ... silly, because there are zero degrees of freedom
##SG> to distinguish "group" from Residual.  It is comforting that the sum of
##SG> the variances sum to the variance of "y", ......
##SG>	 However, I would prefer to have the multilevel software catch this
##SG> case and optionally return an error or drop the redundant group
##SG> with a warning.

## Wrong formula gave a seg.fault at times:
D <-  data.frame(y= rnorm(20,10), ff = gl(4,5),
                 x1=rnorm(20,3), x2=rnorm(20,7))
m0 <- lmer(y ~ (x1 + x2)|ff, data = D)
m1 <- lmer(y ~ x1 + x2|ff  , data = D)
m2 <- lmer(y ~ x1 + (x2|ff), data = D)
m3 <- lmer(y ~ (x2|ff) + x1, data = D)
stopifnot(identical(ranef(m0), ranef(m1)),
          identical(ranef(m2), ranef(m3)),
          inherits(tryCatch(lmer(y ~ x2|ff + x1, data = D), error = function(e)e),
                   "error"))

## Reordering of grouping factors should not change the internal structure
Pm1 <- lmer(strength ~ (1|batch) + (1|sample), Pastes, doFit = FALSE)
Pm2 <- lmer(strength ~ (1|sample) + (1|batch), Pastes, doFit = FALSE)
## The environments of Pm1 and Pm2 should be identical except for
## "call" and "frame":
if(exists("all.equal.X")) {
    all.equal.notCall <- function(x,y,...)
    all.equal.X(env(x), env(y), except = c("call", "frame"), ...)
} else {
    trunclist <- function(x) {
	ll <- as.list(env(x))
	ll[-match(c("call", "frame"), names(ll))]
    }
    all.equal.notCall <- function(x,y,...)
        all.equal(trunclist(x), trunclist(x), ...)
}
stopifnot(all.equal.notCall(Pm1, Pm2))

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
