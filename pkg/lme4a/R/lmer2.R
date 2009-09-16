## Also suppress the warning
assign("det_CHMfactor.warn", TRUE, envir = Matrix:::.MatrixEnv)

##' Install various objects, including Lambda, Lind, Zt and Ut in
##' environment rho based on the list bars of random effects terms and
##' model frame fr.
##'
##' @param bars a list of parsed random-effects terms
##' @param fr a model frame in which to evaluate these terms
##' @rho an environment that is modified by this function
##' @return NULL - the side effect of the function is to modify rho
makeZt <- function(bars, fr, rho)
{
    stopifnot(is.list(bars), all(sapply(bars, is.language)),
              inherits(fr, "data.frame"), is.environment(rho))
    if (!length(bars))
        stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    blist <- lapply(bars,
                    function(x)
                {
                    ff <- eval(substitute(factor(fac), list(fac = x[[3]])), fr)
                    nl <- length(levels(ff))
                    mm <- model.matrix(eval(substitute( ~ foo,
                                                       list(foo = x[[2]]))), fr)
                    nc <- ncol(mm)
                    nseq <- seq_len(nc)
                    sm <- as(ff, "sparseMatrix")
                    if (nc  > 1) 
                        sm <- do.call(rBind, lapply(nseq, function(i) sm))
                    sm@x[] <- t(mm[])
                    list(ff = ff, sm = sm, nc = nc, nl = nl,
                         cnms = colnames(mm))
                })
    nl <- sapply(blist, "[[", "nl")     # no. of levels per term
    ## order terms stably by decreasing number of levels in the factor
    if (any(diff(nl)) > 0) blist <- blist[rev(order(nl))]
    rho$Zt <- do.call(rBind, lapply(blist, "[[", "sm"))
    q <- nrow(rho$Zt)

    ## Create and install Lambda, Lind, etc.  This must be done after
    ## any potential reordering of the terms.
    rho$cnms <- lapply(blist, "[[", "cnms")
    nc <- sapply(blist, "[[", "nc")     # no. of columns per term
    nth <- as.integer((nc * (nc+1))/2)  # no. of parameters per term
    nb <- nc * nl                     # no. of random effects per term
    stopifnot(sum(nb) == q)

    rho$u <- numeric(q)
    rho$theta <- numeric(sum(nth))
    boff <- cumsum(c(0L, nb))           # offsets into b
    thoff <- cumsum(c(0L, nth))         # offsets into theta
    rho$Lambda <-
        do.call(sparseMatrix,
                do.call(rBind,
                        lapply(seq_along(blist), function(i)
                           {
                               mm <- matrix(seq_len(nb[i]), nc = nc[i])
                               dd <- diag(nc[i])
                               ltri <- lower.tri(dd, diag = TRUE)
                               ii <- row(dd)[ltri]
                               jj <- col(dd)[ltri]
                               dd[cbind(ii, jj)] <- seq_along(ii)
                               data.frame(i = as.vector(mm[, ii]) + boff[i],
                                          j = as.vector(mm[, jj]) + boff[i],
                                          x = as.double(rep.int(seq_along(ii),
                                          rep.int(nl[i], length(ii))) +
                                          thoff[i]))
                           })))
    rho$Lind <- as.integer(rho$Lambda@x)
    if (rho$diagonalLambda <- all(nc == 1))
        rho$Lambda <- Diagonal(x = rho$Lambda@x)
    rho$Ut <- crossprod(rho$Lambda, rho$Zt)
    
    lower <- -Inf * (rho$theta + 1)
    lower[unique(diag(rho$Lambda))] <- 0
    rho$lower <- lower
                                        # massage the factor list
    fl <- lapply(blist, "[[", "ff")
    asgn <- seq_along(fl)
                                        # check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        fl <- fl[match(ufn, fnms)]
        asgn <- match(fnms, ufn)
    }
    names(fl) <- ufn
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn
    rho$flist <- fl

    ## Need .f during transition of determinant(L) definition
    rho$.f <- if(package_version(packageDescription("Matrix")$Version) >=
                 "0.999375-31") 2 else 1
    NULL
}

lmer2 <-
    function(formula, data, family = gaussian, REML = TRUE, sparseX = FALSE,
             control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
             compDev = if (sparseX) FALSE else TRUE, subset, weights,
             na.action, offset, contrasts = NULL, ...)
{
    mf <- mc <- match.call()
    if (!missing(family)) {      # call glmer if family is not missing
        mc[[1]] <- "glmer"
        eval(mc, parent.frame())
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)
    if (missing(data)) data <- environment(formula)
                                        # environment for deviance evaluation
    rho <- new.env(parent = environment(formula))
    rho$REML <- as.logical(REML)[1]
                                        # evaluate and install the model frame
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula)
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    rho$frame <- fr <- eval(mf, parent.frame())
    rho$n <- n <- nrow(fr)
                                        # components of the model frame
    y <- model.response(fr)
    if(length(dim(y)) == 1) { # avoid problems with 1D arrays, but keep names
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
    rho$y <- y <- unname(as.double(y))

    rho$weights <- as.numeric(model.weights(fr))
    if (any(rho$weights < 0))
        stop(gettext("negative weights not allowed", domain = "R-lme4"))
    loff <- length(offset <- as.numeric(model.offset(fr)))
    if (loff) {
        if (loff == 1) {
            offset <- rep.int(offset, n)
        } else if (loff != n) {
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                          loff, n), domain = "R-lme4")
        }
    }
    rho$offset <- offset
    
    fe.form <- formula           # evaluate fixed-effects model matrix
    nb <- nobars(formula[[3]])   # fixed-effects terms only
    if (is.null(nb)) nb <- 1
    fe.form[[3]] <- nb
    X <- if (sparseX) {
        sparse.model.matrix(fe.form, fr, contrasts)
    } else as(model.matrix(fe.form, fr, contrasts), "dgeMatrix")
    rownames(X) <- NULL
    rho$X <- X
    rho$sparseX <- sparseX
    
    p <- ncol(X)
    stopifnot((rho$nmp <- n - p) > 0)
    fixef <- numeric(p)
    names(fixef) <- colnames(X)
    rho$fixef <- fixef
    
    attr(fr, "terms") <- NULL
    
    ## Check for method argument which is no longer used
    if (!is.null(method <- list(...)$method)) {
        msg <- paste("Argument", sQuote("method"),
                     "is deprecated.\nUse", sQuote("nAGQ"),
                     "to choose AGQ.  PQL is not available.")
        if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
            warning(msg)
        } else stop(msg)
    }
    
    makeZt(expandSlash(findbars(formula[[3]])), fr, rho)
    
    rho$RX <- chol(rho$XtX <- crossprod(X)) # check for full column rank
    rho$Xty <- unname(as.vector(crossprod(X, y)))
    
    rho$RZX <- rho$Ut %*% X
    rho$fitted <- numeric(n)
    rho$prss <- 0
    rho$ldL2 <- 0
    rho$ldRX2 <- 0
    
    rho$L <- Cholesky(tcrossprod(rho$Zt), LDL = FALSE, Imult = 1)
    ## Zty is stored as a Matrix, not a vector because of a
    ## peculiarity in crossprod(Lambda, Zty) if Lambda is zero
    rho$Zty <- rho$Zt %*% y            
    rho$ZtX <- rho$Zt %*% X
    rho$compDev <- compDev
    sP <- function(x) {
### FIXME: weights are not yet incorporated (needed here?)
        if (compDev) {
            .Call(lme4a:::lmerenv_deviance, parent.env(environment()), x)
        } else {
            stopifnot(length(x) == length(theta))
            theta <<- as.numeric(x)
            Lambda@x[] <<- theta[Lind]           # update Lambda
            Ut <<- crossprod(Lambda, Zt)
            Matrix:::destructive_Chol_update(L, Ut, Imult = 1)
            cu <- solve(L, solve(L, crossprod(Lambda, Zty), sys = "P"),
                        sys = "L")
            RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), sys = "P"),
                          sys = "L")
            RX <<- chol(XtX - crossprod(RZX))
            cb <- solve(t(RX), Xty - crossprod(RZX, cu))
            fixef[] <<- solve(RX, cb)@x
            u[] <<- solve(L, solve(L, cu - RZX %*% fixef, sys = "Lt"),
                          sys = "Pt")@x
            fitted[] <<- (if(length(offset)) offset else 0) +
                (crossprod(Ut, u) + X %*% fixef)@x
            prss <<- sum(c(y - fitted, u)^2) # penalized residual sum of squares
            ldL2[] <<- .f * determinant(L)$mod
            ldRX2[] <<- 2 * determinant(RX)$mod
            if (!REML) return(ldL2 + n * (1 + log(2 * pi * prss/n)))
            ldL2 + ldRX2 + nmp * (1 + log(2 * pi * prss/nmp))
        }
    }
    gP <- function() theta
    gB <- function() cbind(lower = lower,
                           upper = rep.int(Inf, length(theta)))
    environment(sP) <- environment(gP) <- environment(gB) <- rho
    me <- new("lmerenv", setPars = sP, getPars = gP, getBounds = gB)
    if (doFit) {                        # perform the optimization
### FIXME: Allow for selecting an optimizer.  Use optimize for scalar
### problems 
        if (verbose) control$trace <- 1
        nlminb(is.finite(rho$lower), me@setPars, lower = rho$lower,
               control = control)
    }
    me
}

setMethod("fixef", "merenv", function(object, ...) env(object)$fixef)

setMethod("ranef", "merenv", function(object, ...)
          {
              with(env(object), (Lambda %*% u)@x)
          })

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

setMethod("ranef", signature(object = "merenvtrms"),
          function(object, postVar = FALSE, drop = FALSE, whichel = names(ans), ...)
      {
          rho <- env(object)
          ## evaluate the list of matrices
          levs <- lapply(fl <- rho$flist, levels)
          asgn <- attr(fl, "assign")
          nc <- sapply(cnms <- rho$cnms, length)
          nb <- nc * (nl <- unname(sapply(levs, length))[asgn])
          ml <- split(as.vector(rho$Lambda %*% rho$u),
                      rep.int(seq_along(nb), nb))
          for (i in seq_along(ml))
              ml[[i]] <- matrix(ml[[i]], nc = nc[i],
                                dimnames = list(NULL, cnms[[i]]))
          ## produce a list of data frames corresponding to
          ## factors, not terms
          ans <- lapply(seq_along(fl),
                        function(i)
                        data.frame(do.call(cbind, ml[asgn == i]),
                                   row.names = levs[[i]],
                                   check.names = FALSE))
          names(ans) <- names(fl)

          ## Process whichel
          stopifnot(is(whichel, "character"))
          whchL <- names(ans) %in% whichel
          ans <- ans[whchL]

          if (postVar) {
              stop("code not yet written")
### the desired calculation is a diagonal block of
### sigma^2 Lambda(theta)P'L^{-T}L^{-1} P Lambda(theta)
### rewrite this in a general form
          }
          if (drop)
              ans <- lapply(ans, function(el)
                        {
                            if (ncol(el) > 1) return(el)
###                            pv <- drop(attr(el, "postVar"))
                            el <- drop(as.matrix(el))
###                            if (!is.null(pv))
###                                attr(el, "postVar") <- pv
                            el
                        })
          class(ans) <- "ranef.mer"
          ans
      })

devcomp <- function(x, theta, ...)
{
    stopifnot(is(x, "lmerenv"))
    if (!missing(theta)) x@setPars(theta)
    with(env(x),
         list(cmp = c(ldL2 = ldL2, ldRX2 = ldRX2, prss = prss,
              deviance = ldL2 + n * (1 + log(2 * pi * prss/n)),
              REML = ldL2 + ldRX2 + nmp * (1 + log(2 * pi * prss/nmp))),
              dims = c(n = n, p = length(fixef), nmp = nmp, q = nrow(Zt))))
}

setMethod("sigma", signature(object = "lmerenv"),
          function (object, ...) {
              dc <- devcomp(object)
              sqrt(dc$cmp["prss"]/
                   (if (env(object)$REML) dc$dims["nmp"] else dc$dims["n"]))
          })

##' Extract the conditional variance-covariance matrix of the fixed
##' effects 
setMethod("vcov", signature(object = "lmerenv"),
	  function(object, ...)
      {
          en <- env(object)
          rr <- as(sigma(object)^2 * chol2inv(as(en$RX, "matrix")),
                   "dpoMatrix")
          nms <- colnames(en$X)
          dimnames(rr) <- list(nms, nms)
          rr@factors$correlation <- as(rr, "corMatrix")
          rr
      })

setMethod("VarCorr", signature(x = "lmerenv"),
	  function(x, ...)
### Create the VarCorr object of variances and covariances
      {
          sc <- sigma(x)
          rho <- env(x)
          nc <- sapply(cnms <- rho$cnms, length)
          ncseq <- seq_along(nc)
          thl <- split(rho$theta, rep.int(ncseq, (nc * (nc + 1))/2))
          ans <- lapply(ncseq, function(i)
                    {
                        mm <- diag(nrow = nc[i])
                        mm[lower.tri(mm, diag = TRUE)] <- thl[[i]]
                        rownames(mm) <- cnms[[i]]
                        val <- tcrossprod(sc * mm) # variance-covariance
                        stddev <- sqrt(diag(val))
                        correl <- t(val / stddev)/stddev
                        diag(correl) <- 1
                        attr(val, "stddev") <- stddev
                        attr(val, "correlation") <- correl
                        val
                    })
          fl <- rho$flist
          names(ans) <- names(fl)[attr(fl, "assign")]
          attr(ans, "sc") <- sc
          ans
      })

## This is modeled a bit after  print.summary.lm :
printMerenv <- function(x, digits = max(3, getOption("digits") - 3),
                        correlation = TRUE, symbolic.cor = FALSE,
                        signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    dc <- devcomp(x)
    rho <- env(x)
    REML <- rho$REML
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
					       if(getRversion() >= "2.8.0")
					       abbreviate(rn, minlength=1, strict=TRUE)
					       else ## for now
					       .Internal(abbreviate(rn, 1, TRUE))
					       )
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

setMethod("summary", signature(object = "lmerenv"),
	  function(object, ...)
      {
          dc <- devcomp(object)
          rho <- env(object)
          REML <- rho$REML
          fcoef <- fixef(object)
          vcov <- vcov(object)
          corF <- vcov@factors$correlation
          dims <- object@dims
          coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
          llik <- logLik(object, REML)
#          dev <- object@deviance
          mType <- "LMM"
          mName <- switch(mType, LMM = "Linear", NMM = "Nonlinear",
                          GLMM = "Generalized linear",
                          GNMM = "Generalized nonlinear")
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

setMethod("model.frame", signature(formula = "lmerenv"),
	  function(formula, ...) env(formula)$frame)

setMethod("model.matrix", signature(object = "lmerenv"),
	  function(object, ...) env(object)$X)

setMethod("terms", signature(x = "lmerenv"),
	  function(x, ...) attr(env(x)$frame, "terms"))

setMethod("deviance", signature(object="lmerenv"),
	  function(object, REML = NULL, ...)
      {
          if (missing(REML) || is.null(REML) || is.na(REML[1]))
              REML <- env(object)$REML
          devcomp(object)$cmp[if(REML) "REML" else "deviance"]
      })

##' Extract the log-likelihood or restricted log-likelihood
setMethod("logLik", signature(object="lmerenv"),
	  function(object, REML = NULL, ...)
      {
          rho <- env(object)
          if (is.null(REML) || is.na(REML[1]))
              REML <- rho$REML
          val <- -deviance(object, REML = REML)/2
          attr(val, "nall") <- attr(val, "nobs") <- length(rho$y)
          attr(val, "df") <- length(rho$fixef) +
              length(rho$theta) + 1L
          attr(val, "REML") <-  as.logical(REML)
          class(val) <- "logLik"
          val
      })

