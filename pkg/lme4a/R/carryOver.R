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

initial_update <- function(rho, lagged, coef = rep(1, length(lagged[[1]])))
{
    updt <- lapply(lagged, lme4:::accum, coef = coef)
    fl <- rho$flist
    fnms <- names(fl)
    asgn <- attr(fl, "assign")
    stopifnot(length(updt) == 1) # too hard to decide what to do for > 1
    ## for (nm in names(updt)) {
    ##     trm <- which(asgn == match(nm, fnms))
    ##     stopifnot(length(trm) == 1)
        
    ##     rho$trms[[trm]]$A <- rho$trms[[trm]]$Zt <- updt[[nm]]
    ## }
    trm <- which(asgn == match(names(updt), fnms))
    stopifnot(length(trm) == 1)
    trms <- evalbars(rho$formula, rho$frame)
    rows <- diff(rho$rCF@Gp)
    rho
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
    ## Call lmer with doFit = FALSE
    lmerc$idvar <- lmerc$timevar <- lmerc$factors <- NULL
    lmerc$doFit <- FALSE
    lmerc[[1]] <- as.name("lmer")
    lf <- eval.parent(lmerc)
    lst <- as.list(lf)                  # Not sure if we need this
    nolag <- lme4:::merFinalize(lf)
    
    ## create model matrices for the lagged factors
    stopifnot(all(factors %in% names(lf$flist)))
    lagged <- lagged_factor(lf$frame, idvar, timevar, factors)

    return(initial_update(rho, lagged))
   
    ## default is undiscounted model
    dm <- mkZt(updateFL(lf$FL, lagged), nolag@ST)
    undisc <- nolag
    undisc@Zt <- dm$Zt
    undisc@A <- dm$A
    undisc@L <- dm$L
    mer_finalize(undisc)
    undisc@call <- match.call()

    ## fit the AR1 discount formula
    nyp <- length(lagged[[1]]) - 1L     # number of year parameters
    upars <- .Call(mer_ST_getPars, undisc)
    ## parameter bounds
    low <- numeric(length(upars) + 1L)
    up <- c(1, rep(Inf, length(upars)))
    AR1dev <- function(epars) {
        foo <- undisc
        foo@Zt <- mkZt(updateFL(lf$FL, lagged, epars[1]^(0:nyp)), nolag@ST)$Zt
        .Call(mer_ST_setPars, foo, epars[-1])
        .Call(mer_update_dev, foo)
    }
    AR1pars <- nlminb(c(rho = 0.5, upars), AR1dev, control = list(trace = 1),
                      lower = low, upper = up)$par
    AR1 <- undisc
    AR1@Zt <- mkZt(updateFL(lf$FL, lagged, AR1pars[1]^(0:nyp)), nolag@ST)$Zt
    .Call(mer_ST_setPars, AR1, AR1pars[-1])
    .Call(mer_update_dev, AR1)
    .Call(mer_update_ranef, AR1)
    .Call(mer_update_mu, AR1)

    ## fit the general coefficients discount formula
    ## parameter bounds
    ypind <- seq_len(nyp)
    low <- numeric(length(upars) + nyp)
    up <- rep(c(1, Inf), c(nyp, length(upars)))
    AR1dev <- function(epars) {
        foo <- undisc
        foo@Zt <- mkZt(updateFL(lf$FL, lagged, c(1, epars[ypind])), undisc@ST)$Zt
        .Call(mer_ST_setPars, foo, epars[-ypind])
        .Call(mer_update_dev, foo)
    }
    genpars <- nlminb(c(AR1pars[1]^ypind, .Call(mer_ST_getPars, AR1)), AR1dev,
                      control = list(trace = 1), lower = low, upper = up)$par
    gen <- undisc
    gen@Zt <- mkZt(updateFL(lf$FL, lagged, c(1, genpars[ypind])), nolag@ST)$Zt
    .Call(mer_ST_setPars, gen, genpars[-ypind])
    .Call(mer_update_dev, gen)
    .Call(mer_update_ranef, gen)
    .Call(mer_update_mu, gen)
    list(nolag = nolag, undisc = undisc,
         AR1pars = AR1pars, AR1 = AR1,
         genpars = genpars, gen = gen)
}
