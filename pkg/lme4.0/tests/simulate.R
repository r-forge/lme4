require(lme4.0)
source(system.file("test-tools.R", package = "Matrix"))# identical3(),

##' Simple Parametric bootstrap Deviance between two models
##' @title parametric bootstrapped deviance
##' @param m0 two (nested) models fitted to the same data
##' @param m1
##' @param nsim: number of replications
##' @return 2 * logLik( fit(m1), fit(m0) )
##' @author Ben Bolker, Feb.2011;  Martin Maechler, Aug.2011
pboot <- function(m0,m1, nsim = 1, seed=NULL) {
  smat <- simulate(m0, nsim=nsim, seed=seed)
  2 * sapply(smat, function(s) logLik(refit(m1,s)) - logLik(refit(m0,s)))
}

set.seed(54321)

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)

s1 <- simulate(fm1)
stopifnot(is.data.frame(s1), ncol(s1)==1, nrow(s1)==length(fitted(fm1)))  ##  1-column data frame
showProc.time()
s2 <- simulate(fm1,10)
stopifnot(is.data.frame(s2), ncol(s2)==10, nrow(s2)==length(fitted(fm1))) ## 10-column data frame
showProc.time()

r0 <- pboot(fm2,fm1, 10)
showProc.time()
summary(r0)

## binomial (non-Bernoulli)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
	     family = binomial, data = cbpp)
gm0 <- update(gm1, . ~. -period)

s3 <- simulate(gm1)
stopifnot(is.data.frame(s3), ncol(s3[[1]])==2, nrow(s3)==length(fitted(gm1)))
s4 <- simulate(gm1,10)
stopifnot(is.list(s4), sapply(s4,ncol)==2,
	  sapply(s4,nrow) == nrow(cbpp))
showProc.time()


## FIXME: fails on more recent R 3.0.0, maybe 2.15.3?
##  with downdated X'X problem ...
##  17 April 2013; not a problem anymore?
r1 <- pboot(gm0,gm1, 10)
showProc.time()
summary(r1)


## FIXME: want real Poisson example, but will have to simulate one instead for now
nobs <- 50
f <- gl(5,10)
x <- runif(nobs,max=4)
u <- rnorm(5,sd=4)
beta <- c(1,2)
d <- data.frame(f,x)
eta <- model.matrix(~x,data=d) %*% beta + u[f]
mu <- exp(eta)
set.seed(2002)
d$y <- rpois(nobs,lambda=mu)
##

## if uncommented, would require a Suggests: dependency on ggplot2
## library(ggplot2)
##  ggplot(d,aes(x=x,y=y,colour=f))+stat_sum(aes(size=..n..))

gm3 <- glmer(y~x+(1|f),data=d,family=poisson)
s5 <- simulate(gm3,seed=1001)
showProc.time()
stopifnot(is.data.frame(s5), ncol(s5)==1, nrow(s5)==nrow(d))
s6 <- simulate(gm3,10)
showProc.time()
stopifnot(is.data.frame(s6), nrow(s6)==nrow(d), ncol(s6)==10)

invisible(refit(gm3,s6[,1]))

## simulate with offset
d$offset <- rep(1,nobs)
eta <- model.matrix(~x,data=d)%*%beta+u[f]+d$offset
mu <- exp(eta)
set.seed(2002)
d$y <- rpois(nobs,lambda=mu)

gm4 <- glmer(y~x+(1|f),offset=offset,data=d,family=poisson)
s7 <- simulate(gm4,seed=1001)
stopifnot(is.data.frame(s7), nrow(s7)==nrow(d),ncol(s7)==1)
showProc.time()
