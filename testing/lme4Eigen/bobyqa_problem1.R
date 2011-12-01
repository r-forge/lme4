## test file: demonstrates that something is getting changed inside the environment of 'ff'
##  during the course of an unsuccessful bobyqa fit.  When we attempt the same fit
##  a second time we get an error about an inadmissible starting value.

## How do we (1a) make a 'deep copy' of ff (replicating the environment and breaking
## the association with the environment) and (1b) check the contents of the environment
## to see what has changed?  Alternately, can (2) we lock the environment and prevent its
## contents from being changed (hopefully this would lead to some kind of meaningful
## error, if the change in the environment is necessary to the operation, OR would
## allow us to complete the process without modifying/mangling the environment ...)
## Option (3) [ouch] is simply to delve down through levels of the code to see directly
## what's happening ...

library(lme4Eigen)

data(Contraception, package="mlmRev")
Contraception <- within(Contraception, ch <- factor(ifelse(livch != 0, "Y", "N")))

## extract the objective function for the second stage optimization
(ff <- glmer(use ~ age + I(age^2) + ch + (1|district:urban), Contraception, binomial, devFunOnly=2L))
val <- c(get("pp", environment(ff))$theta, get("pp", environment(ff))$beta(0))

## lockEnvironment(environment(ff),bindings=TRUE)  ## attempt to lock environment (doesn't seem to matter)

require(optimx)
## can't run bobyqa twice in a row -- something goes haywire

lo <- c(0, rep(-Inf, 4L)
opt1 <-optimx(val, ff, lower=lo, method="bobyqa", control=list(trace=2L, kkt=FALSE))

## second try
opt1 <-optimx(val, ff, lower=lo, method="bobyqa", control=list(trace=2L, kkt=FALSE))
