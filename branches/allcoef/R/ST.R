checkSTform <- function(ST, STnew)
### Check that the 'STnew' argument matches the form of ST.
{
    stopifnot(is.list(STnew), length(STnew) == length(ST),
              all.equal(names(ST), names(STnew)))
    lapply(seq_along(STnew), function (i)
           stopifnot(class(STnew[[i]]) == class(ST[[i]]),
                     all.equal(dim(STnew[[i]]), dim(ST[[i]]))))
    all(unlist(lapply(STnew, function(m) all(diag(m) > 0))))
}

setMethod("create_A", signature(rCF = "reCovFac"),
          function(rCF, rho, ...)
      {
          rho$A <- t(getLambda(rCF)) %*% rho$Zt
      })

setMethod("create_A", signature(rCF = "ST"),
          function(rCF, rho, ...) .Call(ST_create_A, rCF, rho))

setMethod("getPars", signature(rCF = "ST"),
          function(rCF, ...) .Call(ST_getPars, rCF))

setMethod("getBounds", signature(rCF = "ST"),
          function(rCF, ...) .Call(ST_bounds, rCF))

setMethod("getLambda", signature(rCF = "ST"),
          function(rCF, ...) .Call(ST_Lambda, rCF))

setMethod("setPars", signature(rCF = "ST"),
          function(rCF, pars, rho, ...) .Call(ST_setPars, rCF, pars, rho))


