setMethod("create_A", signature(rCF = "reCovFac"),
          function(rCF, rho, ...) rho$A <- t(getLambda(rCF)) %*% rho$Zt)

setMethod("create_A", signature(rCF = "ST"),
          function(rCF, rho, ...) .Call(ST_create_A, rCF, rho))

setMethod("getPars", signature(x = "ST"),
          function(x, ...) .Call(ST_getPars, x))

setMethod("getBounds", signature(x = "ST"),
          function(x, ...) .Call(ST_bounds, x))

setMethod("getLambda", signature(rCF = "ST"),
          function(rCF, ...) .Call(ST_Lambda, rCF))

setMethod("setPars", signature(x = "ST"),
          function(x, pars, rho, ...) {
              .Call(ST_setPars, x, pars, rho)
              .Call(mer_PIRLS, rho)
          })

setMethod("ranef", signature(object = "ST"),
          function(object, ...)
      {
          dots <- list(...)
          if (is.null(u <- dots$u) || is.null(perm <- dots$perm))
              stop(gettextf("ranef for ST objects requires arguments u and perm"))
          ans <- .Call(ST_create_ranef, object, u, perm)
          ST <- object@ST
          for (i in seq_along(ans))
              colnames(ans[[i]]) <- colnames(ST[[i]])
          ans
      })

setMethod("chol", signature(x = "ST"),
          function(x, ...) .Call(ST_chol, x))

