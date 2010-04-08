eps <- .Machine$double.eps
oneMeps <- 1 - eps
set.seed(1)
etas <-
    lapply(list(-8:8,             # equal spacing to asymptotic area
                runif(20, -8, 8), # random sample from wide uniform dist
                rnorm(20, 0, 8), # random sample from wide normal dist
                -10^30, rnorm(10, 0, 4), 10^30), as.numeric)
etapos <-
    lapply(list(1:20,
                rexp(20),
                rgamma(20, 3),
                pmax(eps, rnorm(20, 2, 1))), as.numeric)

mubinom <-
    lapply(list(runif(100, 0, 1),
                rbeta(100, 1, 3),
                rbeta(100, 0.1, 3),
                c(eps, rbeta(100, 3, 0.1), oneMeps)), as.numeric)

tst.fam <- function(fam, lst, Rname, Cname) 
    lapply(lst, function(x)
           checkEquals(fam[[Rname]](x),
                       .Call(Cname, fam, x, PACKAGE = "lme4a")))
tst.lnki <- function(fam, lst) tst.fam(fam, lst, "linkinv", "family_linkinv")
tst.link <- function(fam, lst) tst.fam(fam, lst, "linkfun", "family_link")
tst.muEta <- function(fam, lst) tst.fam(fam, lst, "mu.eta", "family_muEta")

test.uncons.lnki.R <- function() {  # check link for unconstrained eta
    tst.lnki(binomial(), etas)     # binomial with default, logit link
    tst.muEta(binomial(), etas)
    tst.lnki(binomial("probit"), etas)   # binomial with probit link
    tst.muEta(binomial("probit"), etas)
    tst.lnki(poisson(), etas)   # Poisson with default, log link
    tst.muEta(poisson(), etas)
    tst.lnki(gaussian(), etas)  # Gaussian with default, identity link
    tst.muEta(gaussian(), etas)
}

test.pos.lnki.R <- function() {         # check link for positive eta only
    set.seed(1)
    tst.lnki(Gamma(), etapos)           # gamma family
    tst.muEta(Gamma(), etapos)
    tst.lnki(inverse.gaussian(), etapos) # inverse Gaussian
    tst.muEta(inverse.gaussian(), etapos)    
}

test.binom.link.R <- function() {      # check linkinv for binomial mu
    tst.link(binomial(), mubinom)
    tst.link(binomial("probit"), mubinom)
}

test.pos.link.R <- function() {         # check linkinv for positive mu
    tst.link(poisson(), etapos)
    tst.link(Gamma(), etapos)
    tst.link(inverse.gaussian(), etapos)    
}

test.uncons.link.R <- function() {      # check linkinv for unconstrained mu
    tst.link(gaussian(), etas)
}


