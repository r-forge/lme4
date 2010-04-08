test.link.R <- function() { # check link, inverse link and muEta evaluations
    set.seed(1)
    etas <-
        lapply(list(-8:8,             # equal spacing to asymptotic area
                    runif(20, -8, 8), # random sample from wide uniform dist
                    rnorm(20, 0, 8), # random sample from wide normal dist
                    -10^30, rnorm(10, 0, 4), 10^30), as.numeric)
                 
    tst.lnki <- function(fam) 
        lapply(etas, function(eta)
               checkEquals(fam$linkfun(eta),
                           .Call("family_linkinv", fam, eta, PACKAGE = "lme4a")))
    tst.lnki(binomial())           # binomial with default, logit link
    tst.lnki(binomial("probit"))   # binomial with probit link
    tst.lnki(poisson())            # Poisson with default, log link
    tst.lnki(gaussian())           # Gaussian with default, identity link
## need different selection of etas to test Gamma and inverse.gaussian
}
