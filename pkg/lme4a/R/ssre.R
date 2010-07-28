# Special methods for a random intercepts model.  That is, a model
# with a single, simple, scalar random effect

setClass("sssRe",
         representation(theta = "numeric", L = "numeric",
                        Z = "dgCMatrix", u = "numeric",
                        ff = "factor"),
         validity = function(object) {
             q <- ncol(object@Z)
             n <- nrow(object@Z)
             if (length(object@theta) != 1L)
                 return("theta must have length 1")
             if (length(levels(object@ff)) != q ||
                 length(object@L) != q || length(object@u) != q)
                 return("lengths of u, L and levels(ff) must match ncol(Z)")
         })
         
ssre <- function(formula, data, family, control = list(...),
                 verbose = 0L, nAGQ = 1L, doFit = TRUE, subset,
                 weights, na.action, offset, contrasts = NULL,
                 mustart, etastart, ...)
{
    glmr <- mc <- match.call()
    ## use the glmer function to form an merMod object
    glmr[[1]] <- as.name("glmer")
    doFit.orig <- doFit
    glmr$doFit <- FALSE
    obj <- eval(glmr, parent.frame())
    re <- obj@re
    stopifnot(is(re, "reTrms"),
              length(cnms <- re@cnms) == 1L,
              length(cnms[[1]]) == 1L,
              cnms[[1]][1] == "(Intercept)",
              length(theta <- re@theta) == 1L,
              length(lower <- re@lower) == 1L,
              lower[1] == 0,
              all(re@Lind == 1L),
              TRUE
              )
    new("sssRe", theta = re@theta, u = re@u, L = re@L@x,
        Z = t(re@Zt), ff = re@flist[[1]])
}

