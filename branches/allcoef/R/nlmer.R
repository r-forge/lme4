### To Do: Remove the parameters from the frame before installing it
### in the environment

nlmer <- function(formula, data, start = NULL, verbose = FALSE,
                  nAGQ = 1, doFit = TRUE, subset, weights, na.action,
                  contrasts = NULL, model = TRUE, control = list(), ...)
### Fit a nonlinear mixed-effects model
{
    rho <- lme4:::default_rho()
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
    fe.form <- nlform # modify formula to suppress intercept and add pnames
    fe.form[[3]] <- substitute(0 + foo + bar, 
                               list(foo = parse(text = paste(pnames,
                                                collapse = ' + '))[[1]],
                                    bar = lme4:::nobars(formula[[3]])))
    rho$X <- model.matrix(fe.form, fr)
    rownames(rho$X) <- NULL
    rho$fixef <- numeric(ncol(rho$X))
    names(rho$fixef) <- colnames(rho$X)
    rho$fixef[names(start$nlpars)] <- start$nlpars
    lme4:::lmerFactorList(formula, fr, rho, TRUE, TRUE)
    if (!doFit) return(rho)
#    rho$dims["verb"] <- -1L
    lme4:::merFinalize(rho, control, verbose, mc)
}
