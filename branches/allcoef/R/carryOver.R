## Preliminary version of a function to fit a model with carryover for
## the effect of one grouping factor (e.g. teacher) on another
## grouping factor (e.g. student)

##' Create dgTMatrix objects for lagged factors
##'
##' Given a data frame and the names of variables, idvar and timevar,
##' create a wide data frame from idvar and timevar and any variables
##' named in the factors argument.  Create lagged variables of each of
##' the of the variables named in the factors argument and return them
##' as a list of lists of dgTMatrix objects.
##'
##' @param fr - the original data frame
##' @param idvar - character scalar, the name of the id variable
##' @param timevar - character scalar, the name of the time variable
##' @param factors - a character vector of names of factors to lag over time
##' @return a list of lists of dgTMatrix objects
lagged_factor <- function(fr, idvar = "id", timevar = "Time", factors = "inst")
{
    ## extract names and check arguments for consistency
    fnms <- names(fr)
    stopifnot(length(idvar <- as.character(idvar)) == 1,
              length(timevar <- as.character(timevar)) == 1,
              length(factors <- sapply(factors, as.character)) > 0,
              all(c(idvar, timevar, factors) %in% names(fr)),
              all(unlist(lapply(fr[factors], is.factor))))
    ## ensure timevar is a factor
    fr[[timevar]] <- as.factor(fr[[timevar]])
    ## create the wide data frame
    v.names <- setdiff(fnms, c(idvar, timevar))
    wide <- reshape(fr, direction = "wide", idvar = idvar,
                    timevar = timevar, sep = ".", v.names = v.names)
    nts <- seq_len(nt <- length(tlevs <- levels(fr[[timevar]])))
    varying <- attr(wide, "reshapeWide")$varying
    ## lag the variables named as factors
    for (nm in factors) {
        levs <- levels(fr[[nm]])
        vnms <- paste(nm, tlevs, sep = ".")
        ff <- cbind(wide[, vnms],
                    data.frame(NAs = factor(rep(NA, nrow(wide)),
                               levels = levs))[, rep.int(1L, nt - 1)])
        arr <- array(unlist(ff), c(nrow(wide), nt + 1L, nt))[, nts, nts]
        for (face in nts[-1]) {
            nms <- paste(nm, face, ".", tlevs, sep = '')
            for (j in nts) wide[[nms[j]]] <- factor(arr[,j,face], levels = levs)
            varying <- rbind(varying, nms)
        }
    }
    gennms <- outer(factors, nts[-1], paste, sep = ".")
    ## switch to the long data format
    attr(wide, "reshapeWide")$varying <- varying
    attr(wide, "reshapeWide")$v.names <- c(v.names, as.vector(gennms))
    lagged <- reshape(wide)
    ## match to the original order
    lagged <- lagged[match(interaction(fr[[idvar]], fr[[timevar]], drop = TRUE),
                           interaction(lagged[[idvar]], lagged[[timevar]],
                                       drop = TRUE)), ]
    nms <- cbind(factors, gennms)
    aans <- lapply(seq_along(factors), function(i) {
        ans <- lapply(lagged[nms[i,,drop = TRUE]],
                      function(fac) as(Matrix:::fac2sparse(fac, drop = FALSE),
                                       "dgTMatrix"))
        names(ans) <- tlevs
        ans
    })
    names(aans) <- factors
    aans
}

accum <- function(dgTlst, coef = rep(1, length(dgTlst)))
{
    stopifnot(is.numeric(coef),
              length(coef) == length(dgTlst))
    m1 <- dgTlst[[1]]
    as(new("dgTMatrix",
           i = unlist(lapply(dgTlst, slot, name = "i")),
           j = unlist(lapply(dgTlst, slot, name = "j")),
           Dim = m1@Dim,
           Dimnames = m1@Dimnames,
           x = unlist(lapply(seq_along(coef),
           function(i) coef[i] * dgTlst[[i]]@x))),
       "dgCMatrix")
}

updateFL <- function(FL, lagged, coef = rep(1, length(lagged[[1]])))
{
    updt <- lapply(lagged, accum, coef = coef)
    fl <- FL$fl
    fnms <- names(fl)
    asgn <- attr(fl, "assign")
    for (nm in names(updt)) {
        trm <- which(asgn == match(nm, fnms))
        stopifnot(length(trm) == 1)
        FL$trms[[trm]]$A <- FL$trms[[trm]]$Zt <- updt[[nm]]
    }
    FL
}

carryOver <-
    function(formula, data, family = NULL, REML = TRUE, nAGQ = 1,
             control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, idvar = "id", timevar = "Time",
             factors = "sch", ...)
{
    stopifnot(length(formula <- as.formula(formula)) == 3,
              length(idvar <- as.character(idvar)) == 1,
              length(timevar <- as.character(timevar)) == 1,
              length(factors <- sapply(factors, as.character)) > 0,
              all(c(idvar, timevar, factors) %in% all.vars(formula)))
    lmerc <- mc <- match.call()
    ## Call lmer without Znew and pre and with doFit = FALSE
    lmerc$idvar <- lmerc$timevar <- lmerc$factors <- NULL
    lmerc$doFit <- FALSE
    lmerc[[1]] <- as.name("lmer")
    lf <- eval.parent(lmerc)
    stopifnot(all(factors %in% names(lf$FL$fl)))
    
    ## create model matrices for the lagged factors
    lagged <- lagged_factor(lf$fr$mf, idvar, timevar, factors)
    ## default is undiscounted model
    lf$FL <- updateFL(lf$FL, lagged)
   
#    ans <- do.call(if (!is.null(lf$glmFit))
#                   glmer_finalize else lmer_finalize, lf)
    ans@call <- match.call()
    ans
}
