#### lmer, glmer and nlmer plus methods and utilities

## This environment must have a parent to allow for evaluation of the family initializer
default_rho <- function(parent)
{
    rho <- new.env(parent = parent, hash = TRUE)
    rho$nlmodel <- (~I(x))[[2]]
    rho$muEta <- numeric(0)
    rho$pWt <- numeric(0)
    rho$offset <- numeric(0)
    rho$var <- numeric(0)
    rho$sqrtrWt <- numeric(0)
    rho$ghx <- numeric(0)
    rho$ghw <- numeric(0)
    rho$etaGamma <- array(0, c(0L,1L), list(NULL, "x"))
    rho$nlenv <- new.env(parent = emptyenv())
    rho$deviance <-
	c(ML = NA,	# deviance (for ML estimation)
	  REML = NA,	# REML criterion
	  ldL2 = NA,	# 2 * log(det(L))
	  ldRX2 = NA,	# 2 * log(det(RX))
	  sigmaML = 1,	# ML estimate of common scale
	  sigmaREML = NA,# REML estimate of common scale
	  pwrss = -1,	# penalized, weighted RSS
	  disc = NA,	# discrepancy
	  usqr = -1,	# squared length of u
	  wrss = -1,	# weighted residual sum of squares
	  dev = NA,	# deviance
	  llik = NA,	# log-likelihood
	  NULLdev = 0)	# null deviance
    rho$dims <-
	c(LMM  = 0L,	# not a linear mixed model
	  REML = 0L,	# not REML
	  fTyp = 2L,	# default family is "gaussian"
	  lTyp = 5L,	# default link is "identity"
	  vTyp = 1L,	# default variance function is "constant"
	  nest = 0L,	# not nested
	  useSc= 1L,	# default is to use the scale parameter
	  nAGQ = 1L,	# default is Laplace
	  verb = 0L,	# no verbose output
	  mxit = 300L,	# maximum number of iterations
	  mxfn = 900L,	# max. no. funct. eval.
	  cvg  = 0L)	# no optimization yet attempted
    rho
}

### Utilities for parsing the mixed model formula

#' Return the pairs of expressions that are separated by vertical bars
findbars <- function(term)
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    stopifnot(is.call(term))
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

#' Return the formula omitting the pairs of expressions that are
#' separated by vertical bars
nobars <- function(term)
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
	nb <- nobars(term[[2]])
	if (is.null(nb)) return(NULL)
	term[[2]] <- nb
	return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

#' Substitute the '+' function for the '|' function
subbars <- function(term)
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
	term[[2]] <- subbars(term[[2]])
	return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
	term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

#' Return the list of '/'-separated terms in an expression that
#' contains slashes
slashTerms <- function(x)
{
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

#' from a list of length 2 return recursive interaction terms
makeInteraction <- function(x)
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

#' expand any slashes in the grouping factors returned by findbars
expandSlash <- function(bb)
{
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}

##' Is f1 nested within f2?
##'
##' Does every level of f1 occur in conjunction with exactly one level
##' of f2? The function is based on converting a triplet sparse matrix
##' to a compressed column-oriented form in which the nesting can be
##' quickly evaluated.
##'
##' @param f1 factor 1
##' @param f2 factor 2

##' @return TRUE if factor 1 is nested within factor 2

isNested <- function(f1, f2)
{
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    stopifnot(length(f1) == length(f2))
    sm <- as(new("ngTMatrix",
                 i = as.integer(f2) - 1L,
                 j = as.integer(f1) - 1L,
                 Dim = c(length(levels(f2)),
                 length(levels(f1)))),
             "CsparseMatrix")
    all(diff(sm@p) < 2)
}

##' Create factor list and terms list from r.e. terms.
##'
##' Create the list of model matrices from the random-effects terms in
##' the formula and the model frame.
##'
##' @param formula model formula
##' @param mf model frame
##' @param contrasts
##' @param rmInt logical scalar - should the `(Intercept)` column
##'        be removed before creating Zt
##' @param drop logical scalar indicating if elements with numeric
##'        value 0 should be dropped from the sparse model matrices
##'
evalbars <- function(formula, mf, contrasts, rmInt = FALSE)
{
    ## create factor list for the random effects
    bars <- expandSlash(findbars(formula[[3]]))

### Is this check necessary?  Can we use lmer/nlmer to do a straight glm or nls
### fit?  Probably not but should we be able to do so?
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <-
        lapply(bars,
               function(x)
           {
               ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                                     list(fac = x[[3]])), mf)
               im <- as(ff, "sparseMatrix") # transpose of indicators
               stopifnot(isTRUE(validObject(im, test=TRUE)))
               ## evaluate the model matrix, possibly dropping the intercept
               mm <- model.matrix(eval(substitute(if (rmInt) ~ 0 + expr else ~ expr,
                                                  list(expr = x[[2]]))),
                                  mf, contrasts)
               list(f = ff,
                    Zt = drop0(do.call(rBind,
                    lapply(seq_len(ncol(mm)),
                           function(j) {im@x <- mm[,j]; im}))),
                    ST = `dimnames<-`(diag(nrow = ncol(mm), ncol = ncol(mm)),
                    list(colnames(mm), colnames(mm))))
           })
    ## number of random effects for each term
    nlev <- sapply(fl, function(el) length(levels(el$f)))
    ## order terms by decreasing number of levels in the factor but don't
    ## change the order if this is already true
    if (any(diff(nlev)) > 0) fl <- fl[rev(order(nlev))]

    trms <- lapply(fl, "[", -1)
    names(trms) <- NULL
    fl <- lapply(fl, "[[", "f")
    asgn <- seq_along(fl)
    ## check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        ## check that the lengths of the number of levels coincide
        fl <- fl[match(ufn, fnms)]
        asgn <- match(fnms, ufn)
    }
    names(fl) <- ufn
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn
    list(flist = fl, trms = trms,
         nest =  all(sapply(seq_along(fl)[-1],
         function(i) isNested(fl[[i-1]], fl[[i]]))))
}

##' Install model matrices and the ST object
##'
##' Create the Zt and A sparse matrices from the terms list
##'
##' @param trms the terms list component of the value of evalbars
##'
lmerFactorList <- function(trms, rho)
{
    Ztl <- lapply(trms, `[[`, "Zt")
    Zt <- do.call(rBind, Ztl)
    Zt@Dimnames <- vector("list", 2)
    rho$Zt <- Zt

    ST <- new("ST", ST = lapply(trms, `[[`, "ST"),
              Gp = unname(c(0L, cumsum(sapply(Ztl, nrow)))))
    rho$rCF <- ST
    rho$.active <- c(rho$.active, "rCF")

    rho$A <- create_A(ST, rho)
    n <- length(rho$y)
    p <- length(rho$fixef)
    q <- nrow(Zt)
    rho$eta <- numeric(n)
    rho$mu <- numeric(n)
    rho$resid <- numeric(n)
    rho$perm <- 1:q - 1L                # 0-based identity permutation
    rho$u <- numeric(q)
    rho$RZX <- matrix(0, q, p)
    rho$RX <- matrix(0, p, p)
}

### Control parameters for lmer, glmer and nlmer
lmerControl <- function(trace = getOption("verbose"),
                        iter.max = 300L, eval.max = 900L,
                        algorithm = c("both", "Nelder-Mead", "nlminb"))
    list(iter.max = as.integer(iter.max[1]),
         eval.max = as.integer(eval.max[1]),
	 trace = as.integer(trace),
         algorithm = match.arg(algorithm))

## Old version of lmer
lmer2 <-
    function(formula, data, family = gaussian, REML = TRUE,
             control = list(), start = NULL, verbose = FALSE,
             doFit = TRUE,
             subset, weights, na.action, offset, contrasts = NULL,
             mustart, etastart, ...)
{
    rho <- default_rho(environment(formula))
    rho$deviance["ML"] <- NA            # need to force a copy
    if (missing(data)) data <- environment(formula)
    stopifnot(length(formula <- as.formula(formula)) == 3)
                                        # evaluate and install the model frame
    mf <- mc <- match.call()
    m <- match(c("data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula)
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    n <- nrow(fr)
                                        # components of the model frame
    rho$y <- model.response(fr)
    if(length(dim(rho$y)) == 1) { # avoid problems with 1D arrays, but keep names
        nm <- rownames(rho$y)
        dim(rho$y) <- NULL
        if(!is.null(nm)) names(rho$y) <- nm
    }
    rho$frame <- fr                  # may contain redundant variables
    ## No! needed, e.g. in  (log(Days+1) | . ) ex.:  attr(rho$frame, "terms") <- NULL
    rho$nobs <- nrow(fr)
    rho$weights <- as.numeric(as.vector(model.weights(fr)))
    if (length(rho$weights) == 0)
        rho$weights <- rep.int(1, n)
    if (any(rho$weights < 0))
        stop(gettext("negative weights not allowed", domain = "R-lme4"))
    loff <- length(rho$offset <- as.numeric(as.vector(model.offset(fr))))
    if (loff) {
        if (loff == 1) {
            rho$offset <- rep.int(rho$offset, n)
        } else if (loff != rho$nobs) {
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                          loff, n), domain = "R-lme4")
        }
    }
                                        # starting values expressed as mu or eta
    rho$mustart <- model.extract(mf, "mustart")
    rho$etastart <- model.extract(mf, "etastart")

    fe.form <- formula           # evaluate fixed-effects model matrix
    nb <- nobars(formula[[3]])   # fixed-effects terms only
    if (is.null(nb)) nb <- 1
    fe.form[[3]] <- nb
    rho$X <- model.matrix(fe.form, fr, contrasts)
    rownames(rho$X) <- NULL
    p <- ncol(rho$X)
    rho$start <- numeric(p)             # needed for family$initialize
    rho$fixef <- numeric(p)
    names(rho$fixef) <- colnames(rho$X)
                                        # evaluate and check family
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    ft <- famType(family)
    rho$dims[names(ft)] <- ft
    rho$family <- family
    eval(family$initialize, rho)
                                        # enforce modes on some vectors
    rho$y <- unname(as.double(rho$y))   # must be done after initialize
    rho$mustart <- unname(as.double(rho$mustart))
    rho$etastart <- unname(as.double(rho$etastart))
    if (exists("n", envir = rho))
        rho$n <- as.double(rho$n)
    if (family$family %in% c("binomial", "poisson"))
        rho$dims["useSc"] <- 0L

    largs <- list(...)
    ## Check for method argument which is no longer used
    if (!is.null(method <- largs$method)) {
        msg <- paste("Argument", sQuote("method"),
                     "is deprecated.\nUse", sQuote("nAGQ"),
                     "to choose AGQ.  PQL is not available.")
        if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
            warning(msg)
        } else stop(msg)
        largs <- largs[names(largs) != "method"]
    }
    if(length(largs))
	warning("the following '...' arguments have  *not* been used: ",
		sub("^list", "", deparse(largs, control=NULL)))

    eb <- evalbars(formula, rho$frame, contrasts) # flist, trms, nest
    rho$dims["nest"] <- eb$nest
    rho$flist <- eb$flist
    lmerFactorList(eb$trms, rho)
    if (!(rho$dims["LMM"] <- ft["fTyp"] == 2 && ft["lTyp"] == 5 && ft["vTyp"] == 1)) {
        if (ft["lTyp"] != 5) rho$muEta <- numeric(n)
        if (ft["vTyp"] != 1) rho$var <- numeric(n)
    }
    if (length(rho$var) || length(rho$weights)) rho$sqrtrWt <- numeric(n)
    if (REML && rho$dims["LMM"]) rho$dims["REML"] <- 1L
    rho$u0 <- numeric(length(rho$u))
                                        # evaluate the control argument
    control <- do.call(lmerControl, as.list(control))
    if (!missing(verbose)) control$trace <- as.integer(verbose[1])
    rho$dims["verb"] <- control$trace
    control$trace <- abs(control$trace) # negative values give PIRLS output
    rho$control <- control
    rho$mc <- mc                        # store the matched call

    rho$bds <- getBounds(rho)
    setPars(rho, getPars(rho))          # one evaluation to check structure
    rho$start[] <- rho$fixef            # ensure start is distinct from fixef
    rho$mustart[] <- rho$mu


    if (!doFit) return(rho)
    merFinalize(rho)
}

merFinalize <- function(rho)
{
    if (rho$control$algorithm %in% c("both", "Nelder-Mead")) {
        t0 <- getPars(rho)
        if (length(t0) < 2) {
            res <- optimize(function(x) setPars(rho, x), c(0, 5))
            res$convergence <- 0
            res$par <- res$minimum
        } else {
            res <- optim(getPars(rho), function(x) setPars(rho, x),
                         control = list(maxit = max(1000, rho$control$iter.max),
                         trace = rho$control$trace))
        }
    }
    if (rho$control$algorithm %in% c("both", "nlminb")) {
        control = rho$control
        control$algorithm <- NULL
        res <- nlminb(getPars(rho), function(x) setPars(rho, x),
                      lower = rho$bds[,1], upper = rho$bds[,2],
                      control = control)
    }
    if (res$convergence != 0)
        warning(res$message)
    setPars(rho, res$par)
    nlmodel <- rho$nlmodel
    rho.lst <- as.list(rho)
    rho.lst$nlmodel <- NULL          # it gets quoted in the conversion
    rho.lst <-
        rho.lst[which(names(rho.lst) %in% slotNames(getClass("mer")))]
    rho.lst$Class <- "mer"
    ans <- do.call(new, rho.lst)
    ans@call <- rho$mc
    ans@nlmodel <- nlmodel
    ans
}

famNms <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
            "poisson", "quasibinomial", "quasipoisson", "quasi")
linkNms <- c("logit", "probit", "cauchit", "cloglog", "identity",
             "log", "sqrt", "1/mu^2", "inverse")
varNms <- c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")

famType <- function(family)
{
    if (is.null(family)) return(c(fTyp = 2L, lTyp = 5L, vTyp = 1L))
    if(!is.list(family))
	stop(gettextf("invalid GLM family: %s", format(family),
		      domain = "R-lme4"))
    if (!(fTyp <- match(family$family, famNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    if (!(lTyp <- match(family$link, linkNms, nomatch = 0)))
        stop(gettextf("unknown link: %s",
                      sQuote(family$link), domain = "R-lme4"))
    vNam <- switch(fTyp,
                   "mu(1-mu)",          # binomial
                   "constant",          # gaussian
                   "mu^2",              # Gamma
                   "mu^3",              # inverse.gaussian
                   "mu",                # poisson
                   "mu(1-mu)",          # quasibinomial
                   "mu",                # quasipoisson
                   family$varfun)       # quasi
    if (!(vTyp <- match(vNam, varNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    c(fTyp = fTyp, lTyp = lTyp, vTyp = vTyp)
}

glmer <-
function(formula, data, family = gaussian, start = NULL,
         verbose = FALSE, nAGQ = 1, doFit = TRUE, subset, weights,
         na.action, offset, contrasts = NULL,
         control = list(), mustart, etastart, ...)
{
    mc <- match.call()
    mc[[1]] <- as.name("lmer2")
    eval(mc, parent.frame())
}

#### Extractors specific to mixed-effects models

setMethod("coef", signature(object = "mer"),
	  function(object, ...)
      {
          if (length(list(...)))
              warning(paste('arguments named "',
                            paste(names(list(...)), collapse = ", "),
                                  '" ignored', sep = ''))
          fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
          ref <- ranef(object)
          val <- lapply(ref, function(x) fef[rep(1, nrow(x)),,drop = FALSE])
          for (i in seq(a = val)) {
              refi <- ref[[i]]
              row.names(val[[i]]) <- row.names(refi)
              nmsi <- colnames(refi)
              if (!all(nmsi %in% names(fef)))
                  stop("unable to align random and fixed effects")
              for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
          }
          class(val) <- "coef.mer"
          val
#          new("coef.mer", val)
       })

setAs("mer", "dtCMatrix", function(from)
### Extract the L matrix
      as(from@L, "sparseMatrix"))


##' Extract the fixed effects

setMethod("fixef", signature(object = "mer"),
          function(object, ...)
          object@fixef)

##' Extract the random effects.
##'
##' Extract the conditional modes, which for a linear mixed model are
##' also the conditional means, of the random effects, given the
##' observed responses.  These also depend on the model parameters.
##'
##' @param object an object that inherits from the \code{\linkS4class{mer}} class
##' @param postVar logical scalar - should the posterior variance be returned
##' @param drop logical scalar - drop dimensions of single extent
##' @param whichel - vector of names of factors for which to return results

##' @return a named list of arrays or vectors, aligned to the factor list

setMethod("ranef", signature(object = "mer"),
          function(object, postVar = FALSE, drop = FALSE, whichel = names(ans), ...)
      {
          ## evaluate the list of matrices
          ml <- ranef(object@rCF, u = object@u, perm = object@perm,
                       postVar = postVar, ...)
          ## produce a list of data frames corresponding to factors, not terms
          fl <- object@flist
          asgn <- attr(fl, "assign")
          ans <- lapply(seq_along(fl),
                        function(i)
                            data.frame(do.call(cbind, ml[asgn == i]),
                                       row.names = levels(fl[[i]]),
                                       check.names = FALSE))
          names(ans) <- names(fl)

          ## Process whichel
          stopifnot(is(whichel, "character"))
          whchL <- names(ans) %in% whichel
          ans <- ans[whchL]

          if (postVar) {
### the desired calculation is a diagonal block of
### sigma^2 Lambda(theta)P'L^{-T}L^{-1} P Lambda(theta)
### rewrite this in a general form
              pV <- .Call(ST_postVar, object@rCF, object@L,
                          object@perm, object@flist, whchL)
### need to multiply by sigma^2
              dd <- object@dims
              sc <- 1
              if (dd["useSc"])
                  sc <- object@deviance[if (dd["REML"]) "sigmaREML" else "sigmaML"]
              sc <- sc * sc
              for (i in seq_along(ans))
                  attr(ans[[i]], "postVar") <- pV[[i]] * sc
          }
          if (drop)
              ans <- lapply(ans, function(el)
                        {
                            if (ncol(el) > 1) return(el)
                            pv <- drop(attr(el, "postVar"))
                            el <- drop(as.matrix(el))
                            if (!is.null(pv))
                                attr(el, "postVar") <- pv
                            el
                        })
          class(ans) <- "ranef.mer"
          ans
      })

print.ranef.mer <- function(x, ...) print(unclass(x), ...)
print.coef.mer <- function(x, ...) print(unclass(x), ...)

setMethod("sigma", signature(object = "mer"),
          function (object, ...) {
              dd <- object@dims
              if (!dd["useSc"]) return(1)
              object@deviance[if (dd["REML"]) "sigmaREML"
                              else "sigmaML"]
          })

setMethod("VarCorr", signature(x = "mer"),
	  function(x, ...)
### Create the VarCorr object of variances and covariances
      {
          sc <- sigma(x)
	  ans <- lapply(chol(x@rCF),
                        function(ch) {
                            val <- crossprod(sc * ch) # variance-covariance
                            stddev <- sqrt(diag(val))
                            correl <- t(val / stddev)/stddev
                            diag(correl) <- 1
                            attr(val, "stddev") <- stddev
                            attr(val, "correlation") <- correl
                            val
                        })
          fl <- x@flist
          names(ans) <- names(fl)[attr(fl, "assign")]
          attr(ans, "sc") <- if (x@dims["useSc"]) sc else NA
          ans
      })

#### Methods for standard extractors for fitted models

setMethod("anova", signature(object = "mer"),
	  function(object, ...)
      {
	  mCall <- match.call(expand.dots = TRUE)
	  dots <- list(...)
	  modp <- if (length(dots))
	      sapply(dots, is, "mer") | sapply(dots, is, "lm") else logical(0)
	  if (any(modp)) {		# multiple models - form table
	      opts <- dots[!modp]
	      mods <- c(list(object), dots[modp])
	      names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
				    as.character)
	      mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE),
					attr, "df"))]
	      calls <- lapply(mods, slot, "call")
	      data <- lapply(calls, "[[", "data")
	      if (any(data != data[[1]]))
		  stop("all models must be fit to the same data object")
	      header <- paste("Data:", data[[1]])
	      subset <- lapply(calls, "[[", "subset")
	      if (any(subset != subset[[1]]))
		  stop("all models must use the same subset")
	      if (!is.null(subset[[1]]))
		  header <-
		      c(header, paste("Subset", deparse(subset[[1]]), sep = ": "))
	      llks <- lapply(mods, logLik, REML = FALSE)
	      Df <- sapply(llks, attr, "df")
	      llk <- unlist(llks)
	      chisq <- 2 * pmax(0, c(NA, diff(llk)))
	      dfChisq <- c(NA, diff(Df))
	      val <- data.frame(Df = Df,
				AIC = sapply(llks, AIC),
				BIC = sapply(llks, BIC),
				logLik = llk,
				"Chisq" = chisq,
				"Chi Df" = dfChisq,
				"Pr(>Chisq)" = pchisq(chisq, dfChisq, lower = FALSE),
				row.names = names(mods), check.names = FALSE)
	      class(val) <- c("anova", class(val))
              attr(val, "heading") <-
                  c(header, "Models:",
                    paste(rep(names(mods), times = unlist(lapply(lapply(lapply(calls,
                                           "[[", "formula"), deparse), length))),
                         unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
                         sep = ": "))
	      return(val)
	  }
	  else { ## ------ single model ---------------------
            if (length(object@muEta))
              stop("single argument anova for GLMMs not yet implemented")
            if (length(object@V))
              stop("single argument anova for NLMMs not yet implemented")

            p <- object@dims["p"]
            ss <- object@fixef
#            ss <- (.Call(mer_update_projection, object)[[2]])^2
            names(ss) <- names(object@fixef)
            asgn <- attr(object@X, "assign")
            terms <- terms(object)
            nmeffects <- attr(terms, "term.labels")
            if ("(Intercept)" %in% names(ss))
              nmeffects <- c("(Intercept)", nmeffects)
            ss <- unlist(lapply(split(ss, asgn), sum))
            df <- unlist(lapply(split(asgn,  asgn), length))
            ## dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
            ms <- ss/df
            f <- ms/(sigma(object)^2)
            ## P <- pf(f, df, dfr, lower.tail = FALSE)
            ## table <- data.frame(df, ss, ms, dfr, f, P)
            table <- data.frame(df, ss, ms, f)
            dimnames(table) <-
              list(nmeffects,
                                        #			c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
                   c("Df", "Sum Sq", "Mean Sq", "F value"))
            if ("(Intercept)" %in% nmeffects)
              table <- table[-match("(Intercept)", nmeffects), ]
            attr(table, "heading") <- "Analysis of Variance Table"
            class(table) <- c("anova", "data.frame")
            table
	  }
      })

if (FALSE) {
    setMethod("confint", signature(object = "mer"),
              function(object, parm, level = 0.95, ...)
              .NotYetImplemented()
              )
}

setMethod("deviance", signature(object="mer"),
	  function(object, REML = NULL, ...)
      {
          if (missing(REML) || is.null(REML) || is.na(REML[1]))
              REML <- object@dims["REML"]
          object@deviance[if(REML) "REML" else "ML"]
      })

setMethod("fitted", signature(object = "mer"),
          function(object, ...)
          napredict(attr(object@frame, "na.action"), object@mu))

setMethod("formula", signature(x = "mer"),
	  function(x, ...)
	  x@call$formula
	  )

setMethod("logLik", signature(object="mer"),
	  function(object, REML = NULL, ...)
### Extract the log-likelihood or restricted log-likelihood
      {
          dims <- object@dims
          if (is.null(REML) || is.na(REML[1]))
              REML <- dims["REML"]
          val <- -deviance(object, REML = REML)/2
          attr(val, "nall") <- attr(val, "nobs") <- length(object@y)
          attr(val, "df") <- length(object@fixef) +
              length(getPars(object@rCF)) + as.logical(dims["useSc"])
          attr(val, "REML") <-  as.logical(REML)
          class(val) <- "logLik"
          val
      })
setMethod("predict", signature(object="mer"),
          function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
                    interval = c("none", "confidence", "prediction"),
                    level = 0.95, type = c("response", "terms"),
                    terms = NULL, na.action = na.pass,
                    pred.var = res.var/weights, weights = 1, ...)
          {
            tt <- terms(object)
            if (missing(newdata) || is.null(newdata)) {
              mm <- X <- model.matrix(object)
              mmDone <- TRUE
              offset <- object@offset
            }
            else {
              Terms <- delete.response(tt)
              m <- model.frame(Terms, newdata, na.action = na.action,
                               xlev = object$xlevels)
              if (!is.null(cl <- attr(Terms, "dataClasses")))
                .checkMFClasses(cl, m)
              X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
              offset <- if (!is.null(off.num <- attr(tt, "offset")))
                eval(attr(tt, "variables")[[off.num + 1]], newdata)
              else if (!is.null(object$offset))
                eval(object$call$offset, newdata)
              mmDone <- FALSE
            }
            n <- length(residuals(object))
            predictor <- drop(X %*% fixef(object))
            if (length(offset))
              predictor <- predictor + offset
            return(predictor)

            interval <- match.arg(interval)
            if (interval == "prediction") {
                if (missing(newdata))
                    warning("Predictions on current data refer to _future_ responses\n")
                if (missing(newdata) && missing(weights)) {
                    w <- weights.default(object)
                    if (!is.null(w)) {
                        weights <- w
                        warning("Assuming prediction variance inversely proportional to weights used for fitting\n")
                    }
                }
                if (!missing(newdata) && missing(weights) && !is.null(object$weights) &&
                    missing(pred.var))
                    warning("Assuming constant prediction variance even though model fit is weighted\n")
                if (inherits(weights, "formula")) {
                    if (length(weights) != 2L)
                        stop("'weights' as formula should be one-sided")
                    d <- if (missing(newdata) || is.null(newdata))
                        model.frame(object)
                    else newdata
                    weights <- eval(weights[[2L]], d, environment(weights))
                }
            }
            type <- match.arg(type)
            if (se.fit || interval != "none") {
                res.var <- if (is.null(scale)) {
                    r <- object$residuals
                    w <- object$weights
                    rss <- sum(if (is.null(w)) r^2 else r^2 * w)
                    df <- n - p
                    rss/df
                }
                else scale^2
                if (type != "terms") {
                    if (p > 0) {
                        XRinv <- if (missing(newdata) && is.null(w))
                            qr.Q(object$qr)[, p1, drop = FALSE]
                        else X[, piv] %*% qr.solve(qr.R(object$qr)[p1,
                                                                   p1])
                        ip <- drop(XRinv^2 %*% rep(res.var, p))
                    }
                    else ip <- rep(0, n)
                }
            }
            if (type == "terms") {
                if (!mmDone) {
                    mm <- model.matrix(object)
                    mmDone <- TRUE
                }
                aa <- attr(mm, "assign")
                ll <- attr(tt, "term.labels")
                hasintercept <- attr(tt, "intercept") > 0L
                if (hasintercept)
                    ll <- c("(Intercept)", ll)
                aaa <- factor(aa, labels = ll)
                asgn <- split(order(aa), aaa)
                if (hasintercept) {
                    asgn$"(Intercept)" <- NULL
                    if (!mmDone) {
                        mm <- model.matrix(object)
                        mmDone <- TRUE
                    }
                    avx <- colMeans(mm)
                    termsconst <- sum(avx[piv] * beta[piv])
                }
                nterms <- length(asgn)
                if (nterms > 0) {
                    predictor <- matrix(ncol = nterms, nrow = NROW(X))
                    dimnames(predictor) <- list(rownames(X), names(asgn))
                    if (se.fit || interval != "none") {
                        ip <- matrix(ncol = nterms, nrow = NROW(X))
                        dimnames(ip) <- list(rownames(X), names(asgn))
                        Rinv <- qr.solve(qr.R(object$qr)[p1, p1])
                    }
                    if (hasintercept)
                        X <- sweep(X, 2L, avx, check.margin = FALSE)
                    unpiv <- rep.int(0L, NCOL(X))
                    unpiv[piv] <- p1
                    for (i in seq.int(1L, nterms, length.out = nterms)) {
                        iipiv <- asgn[[i]]
                        ii <- unpiv[iipiv]
                        iipiv[ii == 0L] <- 0L
                        predictor[, i] <- if (any(iipiv > 0L))
                            X[, iipiv, drop = FALSE] %*% beta[iipiv]
                        else 0
                        if (se.fit || interval != "none")
                            ip[, i] <- if (any(iipiv > 0L))
                                as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii,
                                                     , drop = FALSE])^2 %*% rep.int(res.var,
                                                       p)
                            else 0
                    }
                    if (!is.null(terms)) {
                        predictor <- predictor[, terms, drop = FALSE]
                        if (se.fit)
                            ip <- ip[, terms, drop = FALSE]
                    }
                }
                else {
                    predictor <- ip <- matrix(0, n, 0)
                }
                attr(predictor, "constant") <- if (hasintercept)
                    termsconst
                else 0
            }
            if (interval != "none") {
                tfrac <- qt((1 - level)/2, df)
                hwid <- tfrac * switch(interval, confidence = sqrt(ip),
                                       prediction = sqrt(ip + pred.var))
                if (type != "terms") {
                    predictor <- cbind(predictor, predictor + hwid %o%
                                       c(1, -1))
                    colnames(predictor) <- c("fit", "lwr", "upr")
                }
                else {
                    lwr <- predictor + hwid
                    upr <- predictor - hwid
                }
            }
            if (se.fit || interval != "none")
                se <- sqrt(ip)
            if (missing(newdata) && !is.null(na.act <- object$na.action)) {
                predictor <- napredict(na.act, predictor)
                if (se.fit)
                    se <- napredict(na.act, se)
            }
            if (type == "terms" && interval != "none") {
                if (missing(newdata) && !is.null(na.act)) {
                    lwr <- napredict(na.act, lwr)
                    upr <- napredict(na.act, upr)
                }
                list(fit = predictor, se.fit = se, lwr = lwr, upr = upr,
                     df = df, residual.scale = sqrt(res.var))
            }
            else if (se.fit)
                list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
            else predictor
        })



setMethod("residuals", signature(object = "mer"),
	  function(object, ...)
          napredict(attr(object@frame, "na.action"), object@resid))

setMethod("resid", signature(object = "mer"),
	  function(object, ...)
          napredict(attr(object@frame, "na.action"), object@resid))

### FIXME: This method must be rewritten
setMethod("simulate", "mer",
          function(object, nsim = 1, seed = NULL, ...)
      {
	  if(!is.null(seed)) set.seed(seed)
	  if(!exists(".Random.seed", envir = .GlobalEnv))
	      runif(1)		     # initialize the RNG if necessary
          RNGstate <- .Random.seed
          dims <- object@dims
          etasim <- as.vector(object@X %*% fixef(object)) +  # fixed-effect contribution
              sigma(object) * (as(t(object@A) %*%    # random-effects contribution
                               matrix(rnorm(nsim * length(object@u)), nc = nsim),
                                  "matrix")
                               ## residual contribution
                               + matrix(rnorm(nsim * length(object@y)), nc = nsim))
          if (length(object@V) == 0 && length(object@muEta) == 0)
              return(etasim)
          stop("simulate method for GLMMs and NLMMs not yet implemented")
          })

setMethod("summary", signature(object = "mer"),
	  function(object, ...)
      {
          REML <- object@dims["REML"]
          fcoef <- fixef(object)
          vcov <- vcov(object)
          corF <- vcov@factors$correlation
          dims <- object@dims
          coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
          llik <- logLik(object, REML)
          dev <- object@deviance
          mType <- if((non <- as.logical(length(object@etaGamma)))) "NMM" else "LMM"
          if (gen <- as.logical(length(object@muEta)))
              mType <- paste("G", mType, sep = '')
          mName <- switch(mType, LMM = "Linear", NMM = "Nonlinear",
                          GLMM = "Generalized linear",
                          GNMM = "Generalized nonlinear")
	  if(dims["nAGQ"] == 1)
              method <- "the Laplace approximation"
	  else
	      method <- "the adaptive Gaussian Hermite approximation"
          if (mType == "LMM")
              method <- if(REML) "REML" else "maximum likelihood"

          AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
                                 logLik = c(llik),
                                 deviance = dev["ML"],
                                 REMLdev = dev["REML"],
                                 row.names = "")
          if (is.na(AICframe$REMLdev)) AICframe$REMLdev <- NULL
          varcor <- VarCorr(object)
          REmat <- formatVC(varcor)
          if (is.na(attr(varcor, "sc")))
              REmat <- REmat[-nrow(REmat), , drop = FALSE]

          if (nrow(coefs) > 0) {
              if (!dims["useSc"]) {
                  coefs <- coefs[, 1:2, drop = FALSE]
                  stat <- coefs[,1]/coefs[,2]
                  pval <- 2*pnorm(abs(stat), lower = FALSE)
                  coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
              } else {
                  stat <- coefs[,1]/coefs[,2]
                  ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
                  coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
              }
          } ## else : append columns to 0-row matrix ...
          new("summary.mer",
              object,
              methTitle = paste(mName, "mixed model fit by", method),
              logLik = llik,
              ngrps = sapply(object@flist, function(x) length(levels(x))),
              sigma = sigma(object),
              coefs = coefs,
              vcov = vcov,
              REmat = REmat,
              AICtab= AICframe
              )
      })## summary()

setMethod("model.frame", signature(formula = "mer"),
	  function(formula, ...) formula@frame)

setMethod("model.matrix", signature(object = "mer"),
	  function(object, ...) object@X)

setMethod("terms", signature(x = "mer"),
	  function(x, ...) attr(x@frame, "terms"))

setMethod("update", signature(object = "mer"),
	  function(object, formula., ..., evaluate = TRUE)
      {
	  call <- object@call
	  if (is.null(call))
	      stop("need an object with call slot")
	  extras <- match.call(expand.dots = FALSE)$...
	  if (!missing(formula.))
	      call$formula <- update.formula(formula(object), formula.)
	  if (length(extras) > 0) {
	      existing <- !is.na(match(names(extras), names(call)))
	      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	      if (any(!existing)) {
		  call <- c(as.list(call), extras[!existing])
		  call <- as.call(call)
	      }
	  }
	  if (evaluate)
	      eval(call, parent.frame())
	  else call
      })

setMethod("vcov", signature(object = "mer"),
	  function(object, ...)
### Extract the conditional variance-covariance matrix of the fixed effects
      {
          rr <- as(sigma(object)^2 * chol2inv(object@RX),
                   "dpoMatrix")
          nms <- colnames(object@X)
          dimnames(rr) <- list(nms, nms)
          rr@factors$correlation <- as(rr, "corMatrix")
          rr
      })

setMethod("with", signature(data = "mer"),
	  function(data, expr, ...) {
	      dat <- eval(data@call$data)
	      if (!is.null(na.act <- attr(data@frame, "na.action")))
		  dat <- dat[-na.act, ]
	      lst <- c(list(. = data), data@flist, data@frame, dat)
	      eval(substitute(expr), lst[unique(names(lst))])
	  })

### Show and print methods and utilities for them

formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
{
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"), list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(varc, colnames)), "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x, maxlen) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }, maxlen))
	colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
	cbind(reMat, rBind(corr, rep.int("", ncol(corr))))
    } else reMat
}

## This is modeled a bit after  print.summary.lm :
printMer <- function(x, digits = max(3, getOption("digits") - 3),
                     correlation = TRUE, symbolic.cor = FALSE,
                     signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    REML <- so@dims["REML"]
    llik <- so@logLik
    dev <- so@deviance
    dims <- x@dims

    cat(so@methTitle, "\n")
    if (!is.null(x@call$formula))
        cat("Formula:", deparse(x@call$formula),"\n")
    if (!is.null(x@call$data))
        cat("   Data:", deparse(x@call$data), "\n")
    if (!is.null(x@call$subset))
        cat(" Subset:",
            deparse(asOneSidedFormula(x@call$subset)[[2]]),"\n")
    print(so@AICtab, digits = digits)

    cat("Random effects:\n")
    print(so@REmat, quote = FALSE, digits = digits, ...)

    ngrps <- so@ngrps
    n <- length(x@y)
    cat(sprintf("Number of obs: %d, groups: ", n))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    if (is.na(so@sigma))
	cat("\nEstimated scale (compare to 1):",
            sqrt(exp(so@deviance["lr2"])/n), "\n")
    if (nrow(so@coefs) > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(correlation) {
	    corF <- so@vcov@factors$correlation
	    if (!is.null(corF)) {
		p <- ncol(corF)
		if (p > 1) {
		    rn <- rownames(so@coefs)
		    rns <- abbreviate(rn, minlen=11)
		    cat("\nCorrelation of Fixed Effects:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			corf <- as(corF, "matrix")
			dimnames(corf) <- list(rns,
					       abbreviate(rn, minlength=1, strict=TRUE))
			print(symnum(corf))
		    }
		    else {
			corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				       nc = p,
                                       dimnames = list(rns, abbreviate(rn, minlen=6)))
			corf[!lower.tri(corf)] <- ""
			print(corf[-1, -p, drop=FALSE], quote = FALSE)
		    }
		}
	    }
	}
    }
    invisible(x)
}

setMethod("print", "mer", printMer)
setMethod("show", "mer", function(object) printMer(object))

printNlmer <- function(x, digits = max(3, getOption("digits") - 3),
                       correlation = TRUE, symbolic.cor = FALSE,
                       signif.stars = getOption("show.signif.stars"), ...)
### FIXME: Does nlmer need a separate show method?
{
    dims <- x@dims
    cat("Nonlinear mixed model fit by Laplace\n")
    if (!is.null(x@call$formula))
        cat("Formula:", deparse(x@call$formula),"\n")
    if (!is.null(x@call$data))
        cat("   Data:", deparse(x@call$data), "\n")
    if (!is.null(x@call$subset))
        cat(" Subset:",
            deparse(asOneSidedFormula(x@call$subset)[[2]]),"\n")

    cat("Random effects:\n")
    print(formatVC(VarCorr(x)), quote = FALSE,
          digits = max(3, getOption("digits") - 3))

    cat(sprintf("Number of obs: %d, groups: ", length(x@y)))
    ngrps <- sapply(x@flist, function(x) length(levels(x)))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    cat("\nFixed effects:\n")
    print(x@fixef)
    invisible(x)
}

setMethod("refit", signature(object = "mer", newresp = "numeric"),
          function(object, newresp, ...)
      {
          newresp <- as.double(newresp[!is.na(newresp)])
          stopifnot(length(newresp) == length(object@y))
          object@y <- newresp
         # mer_finalize(object)
      })

## cheap (and faster?) version of Matrix::bdiag() or Matrix:::.bdiag() :
BlockDiagonal <- function(lst)
{
    stopifnot(is(lst, "list"))
    lst <- lapply(lapply(lst, as, Class = "generalMatrix"),
                  as, Class = "TsparseMatrix")
    isSquare <- function(x) nrow(x) == ncol(x)
    stopifnot(all(sapply(lst, isSquare)),
              all(sapply(lst, is, class2 = "dMatrix")))
    if ((nl <- length(lst)) == 1) return(lst[[1]])

    offsets <- c(0L, cumsum(sapply(lst, ncol)))
    new("dgTMatrix", Dim = rep.int(offsets[nl + 1], 2),
        i = unlist(lapply(1:nl, function(i) lst[[i]]@i + offsets[i])),
        j = unlist(lapply(1:nl, function(i) lst[[i]]@j + offsets[i])),
        x = unlist(lapply(lst, slot, "x")))
}

### FIXME: This method needs replacing.
setMethod("expand", signature(x = "mer"),
          function(x, sparse = TRUE, ...)
      {
          ind <- seq_along(ST <- x@ST)

          if (!sparse) {
              elexpand <- function(mat)
                  list(T = new("dtrMatrix", uplo = "L", diag = "U",
                       x = as.vector(mat),
                       Dim = dim(mat), Dimnames = dimnames(mat)),
                       S = Diagonal(x = diag(mat)))
              fl <- x@flist
              if (all(attr(fl, "assign") == ind))
                  names(ST) <- names(fl)
              return(lapply(ST, elexpand))
          }

          ## 'Sparse' case :

          nc <- sapply(ST, ncol)
          if(!all(nc >= 1)) stop("some ST entries lack  ncol(.) >= 1")
          nlev <- diff(x@Gp) %/% nc

          Sblock <- function(i) rep(diag(ST[[i]]), each = nlev[i])
          Smat <- Diagonal(x = unname(unlist(lapply(ind, Sblock))))
          if (max(nc) == 1) {
              Tmat <- Matrix:::.diag2tT(Diagonal(ncol(Smat)), uplo="L")
          } else {
              Tblock <- function(i)
              {
                  if (nc[i] == 1) return(Matrix:::.diag2tT(Diagonal(nlev[i]), uplo="L"))
                  STi <- ST[[i]]
                  nci <- nc[i]
                  lt <- lower.tri(STi)
                  offsets <- (1:nlev[i]) - 1L
                  ij <- nlev[i] * (which(lt, arr.ind=TRUE) - 1)
                  new("dtTMatrix", Dim = rep.int(nlev[i] * nci, 2), uplo = "L", diag = "U",
                      x = rep(STi[lt], each = nlev[i]),
                      i = as.integer(outer(offsets, ij[, 1], "+")),
                      j = as.integer(outer(offsets, ij[, 2], "+")))
              }
              Tmat <- as(as(BlockDiagonal(lapply(ind, Tblock)), "triangularMatrix"),
                         "CsparseMatrix")
          }
          list(sigma = sigma(x), P = as(x@L@perm + 1L, "pMatrix"),
               T = as(Tmat, "CsparseMatrix"), S = Smat)
      })

#### Methods for secondary, derived classes

setMethod("deviance", signature(object = "summary.mer"), function(object) object@deviance)
setMethod("logLik", signature(object = "summary.mer"), function(object) object@logLik)
setMethod("vcov", signature(object = "summary.mer"), function(object) object@vcov)
setMethod("summary", signature(object = "summary.mer"), function(object) object)

#### Methods to produce specific plots

plot.coef.mer <- function(x, y, ...)
{
    varying <- unique(do.call("c",
                              lapply(x, function(el)
                                     names(el)[sapply(el,
                                                      function(col)
                                                      any(col != col[1]))])))
    gf <- do.call("rBind", lapply(x, "[", j = varying))
    gf$.grp <- factor(rep(names(x), sapply(x, nrow)))
    switch(min(length(varying), 3),
           qqmath(eval(substitute(~ x | .grp,
                                  list(x = as.name(varying[1])))), gf, ...),
           xyplot(eval(substitute(y ~ x | .grp,
                                  list(y = as.name(varying[1]),
                                       x = as.name(varying[2])))), gf, ...),
           splom(~ gf | .grp, ...))
}

plot.ranef.mer <- function(x, y, ...)
{
    lapply(x, function(x) {
        cn <- lapply(colnames(x), as.name)
        switch(min(ncol(x), 3),
               qqmath(eval(substitute(~ x, list(x = cn[[1]]))), x, ...),
               xyplot(eval(substitute(y ~ x,
                                      list(y = cn[[1]],
                                           x = cn[[2]]))), x, ...),
               splom(~ x, ...))
    })
}

qqmath.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        hw <- 1.96 * se
        list(ylim = range(y - hw, y + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
        panel.grid(h = -1,v = -1)
        panel.abline(h = 0)
        x <- as.numeric(x)
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        ly <- y - 1.96 * se
        uy <- y + 1.96 * se
        panel.segments(x, y - 1.96*se, x, y + 1.96 * se,
                       col = 'black')
        panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(x) {
        if (!is.null(pv <- attr(x, "postVar"))) {
            cols <- 1:(dim(pv)[1])
            se <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
            nr <- nrow(x)
            nc <- ncol(x)
            ord <- unlist(lapply(x, order)) +
                rep((0:(nc - 1)) * nr, each = nr)
            rr <- 1:nr
            ind <- gl(ncol(x), nrow(x), labels = names(x))
            xyplot(unlist(x)[ord] ~
                   rep(qnorm((rr - 0.5)/nr), ncol(x)) | ind[ord],
                   se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
                   scales = list(y = list(relation = "free")),
                   xlab = "Standard normal quantiles",
                   ylab = NULL, aspect = 1, ...)
        } else {
            qqmath(~values|ind, stack(x),
                   scales = list(y = list(relation = "free")),
                   xlab = "Standard normal quantiles",
                   ylab = NULL, ...)
        }
    }
    lapply(x, f)
}

dotplot.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
        if (is.null(se)) return(list())
        x <- as.numeric(x)
        hw <- 1.96 * as.numeric(se[subscripts])
        list(xlim = range(x - hw, x + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16,
                         horizontal = TRUE, col = dot.symbol$col,
                         lty = dot.line$lty, lwd = dot.line$lwd,
                         col.line = dot.line$col, levels.fos = unique(y),
                         groups = NULL, ...)
    {
        x <- as.numeric(x)
        y <- as.numeric(y)
        dot.line <- trellis.par.get("dot.line")
        dot.symbol <- trellis.par.get("dot.symbol")
        sup.symbol <- trellis.par.get("superpose.symbol")
        panel.abline(h = levels.fos, col = col.line, lty = lty, lwd = lwd)
        panel.abline(v = 0, col = col.line, lty = lty, lwd = lwd)
        if (!is.null(se)) {
            se <- as.numeric(se[subscripts])
            panel.segments( x - 1.96 * se, y, x + 1.96 * se, y, col = 'black')
        }
        panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(x, ...) {
        ss <- stack(x)
        ss$ind <- factor(as.character(ss$ind), levels = colnames(x))
        ss$.nn <- rep.int(reorder(factor(rownames(x)), x[[1]]), ncol(x))
        se <- NULL
        if (!is.null(pv <- attr(x, "postVar")))
            se <- unlist(lapply(1:(dim(pv)[1]), function(i) sqrt(pv[i, i, ])))
        dotplot(.nn ~ values | ind, ss, se = se,
                prepanel = prepanel.ci, panel = panel.ci,
                xlab = NULL, ...)
    }
    lapply(x, f, ...)
}

#### Creating and displaying a Markov Chain Monte Carlo sample from
#### the posterior distribution of the parameters

setMethod("mcmcsamp", signature(object = "mer"),
	  function(object, n = 1, verbose = FALSE, saveb = FALSE, ...)
### Generate a Markov chain Monte Carlo sample from the posterior distribution
### of the parameters in a linear mixed model
      {
          ## in lme4: object@fixef <- fixef(object) # force a copy
          n <- max(1, as.integer(n)[1])
          dd <- object@dims
          q <- length(object@u)
          p <- length(object@fixef)
          ranef <- matrix(numeric(0), nrow = q, ncol = 0)
          if (saveb) ranef <- matrix(object@ranef, nrow = q, ncol = n)
          sigma <- matrix(unname(sigma(object)), nrow = 1,
                          ncol = (if (dd["useSc"]) n else 0))
          ff <- object@fixef
          fixef <- matrix(ff, p, n)
          rownames(fixef) <- names(ff)
          ans <- new("merMCMC",
                     Gp = object@Gp,
                     ST = matrix(.Call(mer_ST_getPars, object), length(getPars(object@rCF)), n),
                     call = object@call,
                     dims = object@dims,
                     deviance = rep(unname(object@deviance["ML"]), n),
                     fixef = fixef,
                     nc = sapply(object@ST, nrow),
                     ranef = ranef,
                     sigma = sigma)
          .Call(mer_MCMCsamp, ans, object)
      })

setMethod("HPDinterval", signature(object = "merMCMC"),
          function(object, prob = 0.95, ...)
      {
          nms <- c("fixef", "ST")
          if (length(object@sigma)) nms <- c(nms, "sigma")
          if (length(object@ranef)) nms <- c(nms, "ranef")
          names(nms) <- nms
          lapply(lapply(nms, slot, object = object),
                 HPDinterval, prob = prob)
      })

setMethod("HPDinterval", signature(object = "matrix"),
          function(object, prob = 0.95, ...)
      {
          if (ncol(object) > nrow(object))
              object <- t(object)
          vals <- apply(object, 2, sort)
          if (!is.matrix(vals))
              stop("object must have nsamp > 1")
          nsamp <- nrow(vals)
          npar <- ncol(vals)
          gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
          init <- 1:(nsamp - gap)
          inds <- apply(vals[init + gap, , drop = FALSE] -
                        vals[init, , drop = FALSE], 2, which.min)
          ans <- cbind(vals[cbind(inds, 1:npar)],
                       vals[cbind(inds + gap, 1:npar)])
          dimnames(ans) <- list(colnames(object), c("lower", "upper"))
          attr(ans, "Probability") <- gap/nsamp
          ans
      })

### FIXME: Watch the names of the variance components here
setMethod("VarCorr", signature(x = "merMCMC"),
          function(x, type = c("raw", "varcov", "sdcorr", "logs"), ...)
      {
          if ("raw" == (type <- match.arg(type))) {
              ST <- t(x@ST)
              colnames(ST) <- paste("ST", 1:ncol(ST), sep = '')
              if (length(x@sigma)) return(cbind(ST, sigma = as.vector(x@sigma)))
              return(ST)
          }
          .Call(merMCMC_VarCorr, x, match(type, c("raw", "varcov", "sdcorr", "logs")))
      })

setMethod("as.matrix", signature(x = "merMCMC"),
          function(x, ...)
          cbind(t(x@fixef), VarCorr(x, ...)))

setMethod("as.data.frame", signature(x = "merMCMC"),
          function(x, row.names = NULL, optional = FALSE, ...)
          as.data.frame(as.matrix(x, ...), row.names = row.names, optional = optional, ...))

setAs("merMCMC", "data.frame", function(from) as.data.frame(from))

aslatticeframe <- function(x, ...)
{
    fr <- as.data.frame(x, ...)
    data.frame(dat = unlist(fr),
               par = gl(ncol(fr), nrow(fr), labels = colnames(fr)),
               iter = rep(1:nrow(fr), ncol(fr)))
}

## FIXME: More care should be taken to avoid duplicate argument names
## in the eventual call to lattice functions. Accumulate the arguments
## in a list and use do.call instead of direct calls.

setMethod("xyplot", signature(x = "merMCMC"),
          function(x, data, ...)
      {
          pfr <- aslatticeframe(x, ...)
          xyplot(dat ~ iter|par, pfr,
                 xlab = "Iteration number", ylab = NULL,
                 scales = list(x = list(axs = 'i'),
                 y = list(relation = "free", rot = 0)),
                 type = c("g", "l"),
                 layout = c(1, length(levels(pfr$par))),
                 strip = FALSE, strip.left = TRUE, ...)
      })

setMethod("densityplot", signature(x = "merMCMC"),
          function(x, data, ...)
          densityplot(~ dat | par, aslatticeframe(x, ...),
                      scales = list(relation = 'free'), ...)
          )

setMethod("qqmath", signature(x = "merMCMC"),
          function(x, data, ...)
          qqmath(~ dat | par, aslatticeframe(x, ...),
                 scales = list(y = list(relation = 'free')), ...)
          )


abbrvNms <- function(gnm, cnms)
### Abbreviate names of columns in grouping factors
### gnm - group name
### cnms - column names
{
    ans <- paste(abbreviate(gnm), abbreviate(cnms), sep = '.')
    if (length(cnms) > 1) {
	anms <- lapply(cnms, abbreviate, minlength = 3)
	nmmat <- outer(anms, anms, paste, sep = '.')
	ans <- c(ans, paste(abbreviate(gnm, minlength = 3),
			    nmmat[upper.tri(nmmat)], sep = '.'))
    }
    ans
}

mcmccompnames <- function(ans, object, saveb, trans, glmer, deviance)
### Mangle the names of the columns of the mcmcsamp result ans
### This operation is common to the methods for "lmer" and "glmer"
{
    gnms <- names(object@flist)
    cnms <- lapply(object@ST, colnames)
    ff <- fixef(object)
    colnms <- c(names(ff), if (glmer) character(0) else "sigma^2",
                unlist(lapply(seq_along(gnms),
                              function(i)
                              abbrvNms(gnms[i],cnms[[i]]))))
    if (trans) {
        ## parameter type: 0 => fixed effect, 1 => variance,
        ##		 2 => covariance
        ptyp <- c(integer(length(ff)), if (glmer) integer(0) else 1:1,
                  unlist(lapply(seq_along(gnms),
                                function(i)
                            {
                                k <- length(cnms[[i]])
                                rep(1:2, c(k, (k*(k-1))/2))
                            })))
        colnms[ptyp == 1] <-
            paste("log(", colnms[ptyp == 1], ")", sep = "")
        colnms[ptyp == 2] <-
            paste("atanh(", colnms[ptyp == 2], ")", sep = "")
    }
    if (deviance) colnms <- c(colnms, "deviance")
### FIXME: this will fail for a mer2 object
    if(saveb) {## maybe better colnames, "RE.1","RE.2", ... ?
        .NotYetImplemented()
        rZy <- object@rZy
        colnms <- c(colnms,
                    paste("b", sprintf(paste("%0",
                                             1+floor(log(length(rZy),10)),
                                             "d", sep = ''),
                                       seq_along(rZy)),
                          sep = '.'))
    }
    colnames(ans) <- colnms
    ans
}

#### Odds and ends

## simulestimate <- function(x, FUN, nsim = 1, seed = NULL, control = list())
## {
##     FUN <- match.fun(FUN)
##     stopifnot((nsim <- as.integer(nsim[1])) > 0,
## 	      inherits(x, "lmer"))
##     if (!is.null(seed)) set.seed(seed)
##     ## simulate the linear predictors
##     lpred <- .Call(mer_simulate, x, nsim)
##     sc <- abs(x@devComp[8])
##     ## add fixed-effects contribution and per-observation noise term
##     lpred <- lpred + drop(x@X %*% fixef(x)) + rnorm(prod(dim(lpred)), sd = sc)

##     cv <- do.call(lmerControl, control)
##     Omega <- x@Omega
##     x@wrkres <- x@y <- lpred[,1]
##     .Call(mer_update_ZXy, x)
##     LMEoptimize(x) <- cv
##     template <- FUN(x)
##     if (!is.numeric(template))
##         stop("simulestimate currently only handles functions that return numeric vectors")
##     ans <- matrix(template, nr = nsim, nc = length(template), byrow = TRUE)
##     colnames(ans) <- names(template)
##     for (i in 1:nsim) {
##         x@wrkres <- x@y <- lpred[,i]
##         x@Omega <- Omega
##         .Call(mer_update_ZXy, x)
##         LMEoptimize(x) <- cv
##         foo <- try(FUN(x))
##         ans[i,] <- if (inherits(foo, "try-error")) NA else foo
##     }
##     ans
## }

hatTrace <- function(x)
{
    .NotYetImplemented()
    stopifnot(is(x, "mer"))
}

ST2Omega <- function(ST)
### Temporary function to convert the ST representation of the
### relative variance-covariance matrix returned by lmer into the
### Omega representation required by lmer
{
    if (nrow(ST) == 1) return(as(1/ST^2, "dpoMatrix"))
    dd <- diag(ST)
    T <- as(ST, "dtrMatrix")
    T@diag <- "U"
    crossprod(solve(T)/dd)
}


## setMethod("simulate", signature(object = "mer"),
## 	  function(object, nsim = 1, seed = NULL, ...)
##       {
## 	  if(!exists(".Random.seed", envir = .GlobalEnv))
## 	      runif(1)		     # initialize the RNG if necessary
## 	  if(is.null(seed))
## 	      RNGstate <- .Random.seed
## 	  else {
## 	      R.seed <- .Random.seed
## 	      set.seed(seed)
## 	      RNGstate <- structure(seed, kind = as.list(RNGkind()))
## 	      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
## 	  }

##           stopifnot((nsim <- as.integer(nsim[1])) > 0,
##                     inherits(object, "lmer"))
## 	  ## similate the linear predictors
## 	  lpred <- .Call(mer_simulate, object, nsim)
## 	  sc <- abs(object@devComp[8])

## 	  ## add fixed-effects contribution and per-observation noise term
## 	  lpred <- as.data.frame(lpred + drop(object@X %*% fixef(object)) +
## 				 rnorm(prod(dim(lpred)), sd = sc))
## 	  ## save the seed
## 	  attr(lpred, "seed") <- RNGstate
## 	  lpred
##       })

## We need to define an S4 print method, since using an S3 print
## method fails as soon as you call print() explicitly, e.g. when
## wanting to specify options.

## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

## setMethod("getFixDF", signature(object="mer"),
## 	  function(object, ...) {
## 	      devc <- object@devComp
## 	      rep(as.integer(devc[1]- devc[2]), devc[2])
## 	  })

## simss <- function(fm0, fma, nsim)
## {
##     ysim <- simulate(fm0, nsim)
##     cv <- list(gradient = FALSE, msMaxIter = 200:200,
## 	       msVerbose = 0:0)
##     sapply(ysim, function(yy) {
## 	.Call(mer_update_y, fm0, yy)
## 	LMEoptimize(fm0) <- cv
## 	.Call(mer_update_y, fma, yy)
## 	LMEoptimize(fma) <- cv
## 	exp(c(H0 = fm0@devComp[["logryy2"]],
## 	      Ha = fma@devComp[["logryy2"]]))
##     })
## }

## setMethod("denomDF", "mer",
##           function(x, ...)
##       {
##           mm <- x@X
##           aa <- attr(mm, "assign")
##           tt <- x@terms
##           if (!isNested(x))
##               return(list(coef = as.numeric(rep(NA, length(x@fixef))),
##                           terms = as.numeric(rep(NA,
##                           length(attr(tt, "order"))))))
##           hasintercept <- attr(tt, "intercept") > 0
##           ## check which variables vary within levels of grouping factors
##           vars <- eval(attr(tt, "variables"), x@frame)
##           fl <- x@flist
##           vv <- matrix(0:0, nrow = length(vars), ncol = length(fl),
##                         dimnames = list(NULL, names(fl)))
##           ## replace this loop by C code.
##           for (i in 1:nrow(ans))        # check if variables vary within factors
##               for (j in 1:ncol(ans))
##                   ans[i,j] <- all(tapply(vars[[i]], fl[[j]],
##                                          function(x) length(unique(x)) == 1))
##           ## which terms vary within levels of which grouping factors?
##           tv <- crossprod(attr(tt, "factors"), !ans)
##           ## maximum level at which the term is constant
##           ml <- apply(tv, 1, function(rr) max(0, which(as.logical(rr))))
##           ## unravel assignment applied to terms
##           ll <- attr(tt, "term.labels")
##           if (hasintercept)
##               ll <- c("(Intercept)", ll)
##           aaa <- factor(aa, labels = ll)
##           asgn <- split(order(aa), aaa)
##           nco <- lapply(asgn, length)   # number of coefficients per term
##           nlev <- lapply(fl, function(x) length(levels(x)))
##           if (hasintercept) asgn$"(Intercept)" <- NULL
##           list(ml = ml, nco = nco, nlev = nlev)
##       })

## Utilities for the fitted mer object
slotsz <- function(obj)
    rev(sort(sapply(slotNames(obj), function(s) object.size(slot(obj, s)))))

slotApply <- function(object, f, ..., simplify = FALSE) {
   .localFun <- function(what, ...) f(slot(object, what), ...)
   sapply(slotNames(object), .localFun, ..., simplify = simplify)
}


yfrm <- function(fm)
{
    stopifnot(is(fm, "mer"))
    snr <- slotApply(fm, function(x)
                 {
                     if (is(x, "matrix") ||
                         is(x, "data.frame") ||
                         is(x, "numeric")) return (NROW(x))
                     0
                 }, simplify = TRUE)
    snr <- snr[snr > 0 & !(names(snr) %in%
                           c("Gp", "dims", "deviance", "frame", "flist", "X"))]
    fr <- cbind(fm@frame, fm@flist[1:NROW(fm@frame), !(names(fm@flist) %in%
                                     names(fm@frame))])
    n <- NROW(fr)
    if (NROW(fm@X) == n)
        fr <- cbind(fr, X = fm@X, Xbeta = fm@X %*% fm@fixef,
                    Zb = crossprod(fm@Zt, fm@ranef)@x)
    do.call(cbind, c(list(fr), sapply(names(which(snr == NROW(fr))),
                                      slot, object = fm, simplify = FALSE)))
}

##' Evaluate conditional components of an LMM.
##'
##' Evaluate conditional components of a linear mixed model for a grid of ST
##' parameter values.
##'
##' @param fm - a fitted linear mixed model
##' @param parmat - a numeric matrix whose rows constitute suitable parameter
##'     values for fm@ST
##' @param type - which slot to extract
##' @return a data frame of deviance values or fixed-effects or random effects
##' @keywords models
devmat <-
    function(fm, parmat, slotname = c("deviance", "fixef", "ranef", "u"), ...)
{
    stopifnot(is(fm, "mer"))
    dd <- fm@dims
    stopifnot(dd["fTyp"] == 2L, # gaussian family
              dd["lTyp"] == 5L, # identity link
              dd["vTyp"] == 1L, # variance function is "constant"
              length(fm@V) == 0L, # nonlinear parameter gradient is identity
              length(fm@muEta) == 0L) # eta -> mu map is identity
    oldpars <- getPars(fm)

    parmat <- as.matrix(parmat)
    storage.mode(parmat) <- "double"
    if (ncol(parmat) == dd["np"])
        parmat <- t(parmat)             # parameter vectors as columns
    stopifnot(nrow(parmat) == dd["np"])
    slotname <- match.arg(slotname)

    slotval <- function(x) {            # function to apply
        setPars(fm, x)
        if (slotname != "deviance") stop("Code not yet written")
        slot(fm, slotname)
    }
    ans <- apply(parmat, 2, slotval)
    slotval(oldpars)                    # restore the fitted model
    as.data.frame(t(rbind(parmat, ans)))
}

if (FALSE) {
### FIXME: Move this function to the stats package
rWishart <- function(n, df, invScal)
### Random sample from a Wishart distribution
    .Call(lme4_rWishart, n, df, invScal)
}

simGLMM <- function(formula, data, family, theta,
                    fixef = rep.int(1, p), verbose = FALSE,
                    control = list(), ...)
{
    rho <- default_rho(environment(formula))
    stopifnot(inherits(data, "data.frame"),
              length(formula <- as.formula(formula)) == 3,
              is.name(ynm <- formula[[2]]))
    ynm <- as.character(ynm)
    if (is.null(data[[ynm]])) data[[ynm]] <- 1

    mf <- mc <- match.call()    # evaluate and install the model frame
    m <- match(c("data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- substitute(~ foo, list(foo = subbars(formula)[[3]]))
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    n <- nrow(fr)

    rho$frame <- fr                  # may contain redundant variables
    attr(rho$frame, "terms") <- NULL
    rho$weights <- as.numeric(as.vector(model.weights(fr)))
    if (length(rho$weights) == 0)
        rho$weights <- rep.int(1, n)
    if (any(rho$weights < 0))
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
    rho$mustart <- model.extract(mf, "mustart")
    rho$etastart <- model.extract(mf, "etastart")

    fe.form <- formula           # evaluate fixed-effects model matrix
    nb <- nobars(formula[[3]])   # fixed-effects terms only
    if (is.null(nb)) nb <- 1
    fe.form <- eval(substitute(~ foo, list(foo = nb)))
    rho$X <- model.matrix(fe.form, fr, contrasts)
    rownames(rho$X) <- NULL
    p <- ncol(rho$X)
    stopifnot(is.numeric(fixef), length(fixef) == p)

    rho$start <- fixef                  # needed for family$initialize
    rho$fixef <- fixef
    names(rho$fixef) <- colnames(rho$X)
    #lmerFactorList(formula, rho$frame, rho, contrasts)
    if (!missing(theta)) {
        theta <- as.double(theta)
        stopifnot(length(theta) == length(theta0 <- getPars(rho$rCF)))
        setPars(rho, theta)
    }
    q <- nrow(rho$Zt)
    rho$u <- rnorm(q)
    rho$eta <- numeric(n)
    rho$y <- numeric(n)
    rho$mu <- numeric(n)
    rho$resid <- numeric(n)
    rho$var <- numeric(n)
    rho$muEta <- numeric(n)
    rho$sqrtrWt <- numeric(n)
    rho$u0 <- numeric(q)
                                       # evaluate and check family
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    ft <- famType(family)
    rho$dims[names(ft)] <- ft
    rho$family <- family
    rho$nobs <- n
    eval(family$initialize, rho)
                                        # enforce modes on some vectors
    .Call("mer_update_mu", rho)
    rho$y <- switch(family$family,
                    poisson = rpois(n, rho$mu),
                    binomial = rbinom(n, 1, rho$mu))
    .Call("mer_update_mu", rho)

    if (!missing(verbose)) control$trace <- as.integer(verbose[1])
    rho$dims["verb"] <- control$trace
    control$trace <- abs(control$trace) # negative values give PIRLS output
    rho$control <- control
    rho$mc <- match.call()

    rho
}

setMethod("evalDev", signature(object = "mer", pars = "matrix"),
          function(object, pars, ...)
      {
          cc <- object@call
          cc$doFit <- FALSE
          cc$verbose <- FALSE
          rho <- eval(cc, parent.frame())
          optpars <- getPars(rho)
          stopifnot(ncol(pars) == length(optpars))
          storage.mode(pars) <- "double"
          cbind(pars, apply(pars, 1, setPars, x = rho))
      })
