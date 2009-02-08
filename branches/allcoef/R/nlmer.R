##' Fit a nonlinear mixed-effects model
##'
##' @param formula a nonlinear mixed model formula (see detailed documentation)
##' @param data an optional data frame containing the variables named in
##'    \code{formula}.  By default the variables are taken from the
##'    environment from which \code{nlmer} is called.
##' @param start starting estimates for the nonlinear model
##'    parameters, as a named numeric vector 
##' @param verbose integer scalar passed to nlminb.  If negative then
##'    diagnostic output from the PIRLS (penalized iteratively
##'    reweighted least squares) step is also provided.
##' @param nAGQ number of adaptive Gauss-Hermite quadrature points to use
##' @param doFit logical scalar.  If FALSE the optimization
##'    environment is returned. Otherwise the parameters are estimated
##'    and an object of S4 class "mer" is returned.
##' @param subset further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param weights  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param na.action  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param contrasts  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param control a list of control parameters passed to nlminb.  The
##'    defaults are given in the (hidden) function \code{lmerControl}.

##' @return if doFit is FALSE an environment, otherwise an object of S4 class "mer"
nlmer <- function(formula, data, start = NULL, verbose = FALSE,
                  nAGQ = 1, doFit = TRUE, subset, weights, na.action,
                  contrasts = NULL, control = list(), ...)
{
    rho <- default_rho(environment(formula))
    mf <- mc <- match.call()
    m <- match(c("data", "subset", "weights", "na.action", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")

    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a 3-part formula")
    nlform <- as.formula(formula[[2]])
    if (length(nlform) < 3)
        stop("formula must be a 3-part formula")
    nlmod <- as.call(nlform[[3]])
                                        # check for parameter names in start 
    if (is.numeric(start)) start <- list(nlpars = start)
    stopifnot((s <- length(pnames <- names(start$nlpars))) > 0,
              is.numeric(start$nlpars))
    if (!all(pnames %in% (anms <- all.vars(nlmod))))
        stop("not all parameter names are used in the nonlinear model expression")
    fr.form <- nlform
    fr.form[[3]] <-         # the frame formula includes all variables
        parse(text = paste(setdiff(all.vars(formula), pnames),
                         collapse = ' + '))[[1]]
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    rho$y <- model.response(fr, "double")
    rho$weights <- as.numeric(as.vector(model.weights(fr)))
    if (length(rho$weights) && any(rho$weights) < 0)
        stop(gettext("negative weights not allowed", domain = "R-lme4"))
    loff <- length(rho$offset <- as.numeric(as.vector(model.offset(fr))))
    if (loff) {
        if (loff == 1) {
            rho$offset <- rep.int(rho$offset, rho$nobs)
        } else if (loff != rho$nobs) {
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                          loff, rho$nobs), domain = "R-lme4")
        }
    }
    rho$frame <- fr
    attr(rho$frame, "terms") <- NULL
    for (nm in pnames) fr[[nm]] <- start$nlpars[[nm]]
    n <- nrow(fr)

    ## create nlenv and check the evaluation of the nonlinear model function
    rho$nlenv <- new.env()  # want it to inherit from this environment (or formula env)
    lapply(all.vars(nlmod), function(nm) assign(nm, fr[[nm]], envir = rho$nlenv))
    rho$nlmodel <- nlmod
    if (is.null(rho$etaGamma <- attr(eval(rho$nlmodel, rho$nlenv), "gradient")))
        stop("The nonlinear model in nlmer must return a gradient attribute")

    ## build the extended frame for evaluation of X and Zt
    fr <- do.call(rbind, lapply(1:s, function(i) fr)) # rbind s copies of the frame
    for (nm in pnames) # convert these variables in fr to indicators
        fr[[nm]] <- as.numeric(rep(nm == pnames, each = n))
    fe.form <- nlform # modify formula to suppress intercept (Is this a good idea?)
    fe.form[[3]] <- substitute(0 + bar, list(bar = nobars(formula[[3]])))
    rho$X <- model.matrix(fe.form, fr)
    rownames(rho$X) <- NULL
    if ((qrX <- qr(rho$X))$rank < ncol(rho$X))
        stop(gettextf("rank of X = %d < ncol(X) = %d", qrX$rank, ncol(rho$X)))
    rho$start <- rho$fixef <- qr.coef(qrX, unlist(lapply(pnames, get, envir = rho$nlenv)))
    lmerFactorList(formula, fr, rho, TRUE, TRUE)
    if (!doFit) return(rho)
    merFinalize(rho, control, verbose, mc)
}
