###################################################
### chunk number 1: preliminaries
###################################################
options(width=69,show.signif.stars=FALSE)
library(lattice)
lattice.options(default.theme = function() standard.theme())
#lattice.options(default.theme = function() standard.theme(color=FALSE))
library(lme4a)
fm1 <- lmer(Yield ~ 1 + (1|Batch), Dyestuff)
fm1ML <- update(fm1, REML = FALSE)
VCREML <- VarCorr(fm1)
VCML <- VarCorr(fm1ML)
##' Return an interpolated minimum on a grid
##'
##' @param x numeric values of dimension being minimized.  Default is seq_along(y)
##' @param y numeric response vector
##' @return location (the x value) of the minimum and its value as
##'         obtained by quadratic interpolation at three values
##'         bracketing the minimum. 
gmin <- function(x = seq_along(y), y)
{
    conv <- array(c(0,-1,1,2,0,-2,0,1,1), c(3,3))/2
    mpos <- which.min(y)
    if (mpos %in% c(1L, length(y))) return(c(x[mpos], y[mpos]))
    ycoef <- conv %*% y[mpos + -1:1]
    pos <- -ycoef[2]/(2*ycoef[3])
    stopifnot(-1 < pos, pos < 1)
    b <- c(1, pos, pos^2)
    c(crossprod(conv %*% x[mpos + -1:1], b), crossprod(ycoef, b))
}
##' Deviance for a specified sigma (that is, without profiling on sigma)
##'
##' @param fenv  evaluation environment
##' @param theta relative covariance parameters
##' @param sigma common scale parameter
##' @return vector of the deviance and the REML criterion
sigmaDev <- function(fenv, theta, sigma)
{
    setPars(fenv, theta)
    sigmasq <- sigma^2
    dev <- fenv$deviance
    n <- nrow(fenv$X)
    nmp <- n - ncol(fenv$X)
    base <- unname(dev["ldL2"] +  dev["pwrss"]/sigmasq)
    l2ps <- log(2*pi*sigmasq)
    base + unname(c(dev["ldRX2"],0)) + c(nmp,n) * l2ps
}
##' Deviance at a specified sigma, alternative argument list
##'
##' @param fenv  evaluation environment
##' @param pars  concatenation of relative covariance parameters and sigma
##' @return vector of the deviance and the REML criterion
##' @comment FIXME, change this to be an option of sigmaDev
sigmaDev1 <- function(fenv, pars)
{
    pars <- as.numeric(pars)
    stopifnot((lp <- length(pars)) > 1)
    sigmaDev(fenv, pars[-lp]/pars[lp], pars[lp])
}

if (file.exists("grd1.rda")) {
    load("grd1.rda")
    sigB <- sort(unique(grd1$sigB))
    sigmas <- sort(unique(grd1$sigma))
    deviance <- array(grd1$deviance + deviance(fm1ML),
                      c(length(sigB), length(sigmas)))
    REML <- array(grd1$REML + deviance(fm1),
                  c(length(sigB), length(sigmas)))
} else {
    fm1env <- lmer(Yield ~ 1 + (1|Batch), Dyestuff, doFit = FALSE)
    n <- nrow(fm1env$X)
    nmp <- n - ncol(fm1env$X)
    sigB <- seq(0, 225,len = 101)
    sigmas <- seq(30, 100,len = 101)
    grd1 <- expand.grid(sigB = sigB, sigma = sigmas)
    attr(grd1, "out.attrs") <- NULL          # save a bit of space
    vals <- apply(t(grd1[, 1:2]), 2, sigmaDev1, fenv = fm1env)
    grd1$REML <- vals[1,] - deviance(fm1)
    grd1$deviance <- vals[2,] - deviance(fm1ML)
    deviance <- array(vals[2,], c(length(sigB), length(sigmas)))
    REML <- array(vals[1,], c(length(sigB), length(sigmas)))
    save(grd1, file = "grd1.rda")
}
sigmasq <- sigmas * sigmas
sigBsq <- sigB * sigB
parvec <- rep.int(c("sigma", "sigmaB"), c(length(sigmas), length(sigB)))
profsd <- t(apply(deviance, 1, gmin, x = sigmas))
profsbd <- t(apply(deviance, 2, gmin, x = sigB))
profsR <- t(apply(REML, 1, gmin, x = sigmas))
profsbR <- t(apply(REML, 2, gmin, x = sigB))
sqrts <- sqrt(profsbd[,2] - deviance(fm1ML))
sqrtsb <- sqrt(profsd[,2] - deviance(fm1ML))
ssqrtb <- sqrtsb * ifelse(seq_along(sqrtsb) < which.min(sqrtsb), -1, 1)
ssqrt <- sqrts * ifelse(seq_along(sqrts) < which.min(sqrts), -1, 1)


###################################################
### chunk number 2: DyeML
###################################################
fm1env <- lmer(Yield ~ 1 + (1|Batch), Dyestuff, REML=FALSE, doFit = FALSE)
n <- length(fm1env$y)
theta <- seq(0,3,0.01)
res <- cbind(as.data.frame(t(sapply(theta,
                                  function(x) {
                                      setPars(fm1env, x)
                                      fm1env$deviance[c(1:7,9,10)]
                                  }))), theta)
names(res)[1] <- "deviance"
res$lprss <- n * (1 + log(2 * pi * res$pwrss/n))
print(xyplot(deviance + ldL2 + lprss ~ theta,
             res, xlab = expression(theta), ylab = NULL, outer = TRUE,
             scales = list(y = list(relation = "free", rot = 0),
                           x = list(axs = 'i')),
             aspect = 1, type = c("g","l"), layout = c(3,1)))


###################################################
### chunk number 3: Dyestuff2dot
###################################################
print(dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff2,
              ylab = "Batch", jitter.y = TRUE, pch = 21,
              xlab = "Simulated response (dimensionless)",
              type = c("p", "a")))


###################################################
### chunk number 4: Dye2ML
###################################################
fm2MLpart <- lmer(Yield ~ 1 + (1|Batch), Dyestuff2, REML=FALSE, doFit = FALSE)
theta <- seq(0,3,0.01)
res2 <- cbind(as.data.frame(t(sapply(theta,
                                  function(x) {
                                      setPars(fm2MLpart, x)
                                      fm2MLpart$deviance[c(1:7,9,10)]
                                  }))), theta)
names(res2)[1] <- "deviance"
res2$lprss <- n * (1 + log(2 * pi * res2$pwrss/n))
print(xyplot(deviance + ldL2 + lprss ~ theta,
             res2, xlab = expression(theta), ylab = NULL, outer = TRUE,
             scales = list(y = list(relation = "free", rot = 0),
                           x = list(axs = 'i')),
	     aspect = 1, type = c("g","l"), layout = c(3,1)))
    


###################################################
### chunk number 5: devcontoursfm1
###################################################
print(contourplot(deviance + REML ~ sigB * sigma, grd1,
                  xlab = expression(sigma[B]),
                  ylab = expression(sigma),
                  aspect = 1,
                  at = 
                  qchisq(c(0.5,0.8,0.9,0.95,0.99,0.999), df = 2),
                  labels = list(labels = paste(c(50,80,90,95,99,99.9),
                                "%", sep = "")),
                  panel = function(...)
              {
                  if (panel.number() == 1) {
                      VC <- VCML
                      prs <- profsd
                      prsb <- profsbd
                  } else {
                      VC <- VCREML
                      prs <- profsR
                      prsb <- profsbR
                  }
                  panel.points(attr(VC[[1]],"stddev"),
                               attr(VC,"sc"), pch = 3)
                  panel.grid(h = -1, v = -1)
                  panel.lines(sigB, prs[,1], lty = 2)
                  panel.lines(prsb[,1], sigmas, lty = 3)                  
                  panel.contourplot(...)
              }
                  ))


###################################################
### chunk number 6: profileddevfm1
###################################################
print(xyplot(dev ~ x|par,
             data.frame(dev = c(profsbd[,2], profsd[,2]),
                        x = c(sigmas, sigB),
                        par = parvec),
             strip = function(..., factor.levels)
             strip.default(..., factor.levels = c(expression(sigma),
                                expression(sigma[B]))),
             ylab = "Deviance", xlab = NULL, type = c("g","l"),
             scales = list(x = list(relation = "free", axs = 'i'))))


###################################################
### chunk number 7: profileddevsqfm1
###################################################
print(xyplot(dev ~ x|par,
             data.frame(dev = c(profsbd[,2], profsd[,2]),
                        x = c(sigmas^2, sigB^2),
                        par = gl(2, length(sigmas))),
             strip = function(..., factor.levels)
             strip.default(..., factor.levels = c(expression(sigma^2),
                                expression(sigma[B]^2))),
             ylab = "Deviance", xlab = NULL, type = c("g","l"),
             scales = list(x = list(relation = "free", axs = 'i'))))


###################################################
### chunk number 8: profzfm1
###################################################
print(xyplot(profz ~ x|par,
             data.frame(profz = c(sqrts, sqrtsb),
                        x = c(sigmasq, sigBsq),
                        par = parvec),
             strip = function(..., factor.levels)
             strip.default(..., factor.levels = c(expression(sigma^2),
                                expression(sigma[B]^2))),
             ylab = "Profile z", xlab = NULL, type = c("g","l"),
             scales = list(x = list(relation = "free", axs = 'i'))))


###################################################
### chunk number 9: sprofzfm1
###################################################
print(xyplot(profz ~ x|par,
             data.frame(profz = c(ssqrt, ssqrtb),
                        x = c(sigmasq, sigBsq),
                        par = parvec),
             strip = function(..., factor.levels)
             strip.default(..., factor.levels = c(expression(sigma^2),
                                expression(sigma[B]^2))),
             ylab = "Profile z", xlab = NULL, type = c("g","l"),
             scales = list(x = list(relation = "free", axs = 'i'))))


