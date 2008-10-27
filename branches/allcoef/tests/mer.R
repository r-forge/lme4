library(lme4)

VecFromNames <- function(nms, mode = "numeric", defaults = list())
### Generate a named vector of the given mode
{
    ans <- vector(mode = mode, length = length(nms))
    names(ans) <- nms
    ans[] <- NA
    if ((nd <- length(defaults <- as.list(defaults))) > 0) {
        if (length(dnms <- names(defaults)) < nd)
            stop("defaults must be a named list")
        stopifnot(all(dnms %in% nms))
        ans[dnms] <- as(unlist(defaults), mode)
    }
    ans
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

##' dimsNames and devNames are in the package's namespace rather than
##' in the function lmerFactorList because the function sparseRasch
##' needs to access them. 

dimsNames <- c("nt", "n", "p", "q", "s", "np", "LMM", "REML",
               "fTyp", "lTyp", "vTyp", "nest", "useSc", "nAGQ",
               "verb", "mxit", "mxfn", "cvg")
dimsDefault <- list(s = 1L,             # identity mechanistic model
                    mxit= 300L,         # maximum number of iterations
                    mxfn= 900L, # maximum number of function evaluations
                    verb= 0L,           # no verbose output
                    np= 0L,             # number of parameters in ST
                    LMM= 0L,            # not a linear mixed model
                    REML= 0L,         # glmer and nlmer don't use REML
                    fTyp= 2L,           # default family is "gaussian"
                    lTyp= 5L,           # default link is "identity"
                    vTyp= 1L, # default variance function is "constant"
                    useSc= 1L, # default is to use the scale parameter
                    nAGQ= 1L,                  # default is Laplace
                    cvg = 0L)                  # no optimization yet attempted
                    
devNames <- c("ML", "REML", "ldL2", "ldRX2", "sigmaML",
              "sigmaREML", "pwrss", "disc", "usqr", "wrss",
              "dev", "llik", "NULLdev")


##' Create model matrices from r.e. terms.
##'
##' Create the list of model matrices from the random-effects terms in
##' the formula and the model frame.
##' 
##' @param formula model formula
##' @param mf model frame
##' @param rmInt logical scalar - should the `(Intercept)` column
##'        be removed before creating Zt
##' @param drop logical scalar indicating if elements with numeric
##'        value 0 should be dropped from the sparse model matrices 
##'
##' @return a list with components named \code{"trms"}, \code{"fl"}
##'        and \code{"dims"}
lmerFactorList <- function(formula, fr, rmInt, drop)
{
    mf <- fr$mf
    ## record dimensions and algorithm settings

    ## create factor list for the random effects
    bars <- lme4:::expandSlash(lme4:::findbars(formula[[3]]))
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars,
                 function(x)
             {
                 ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                                       list(fac = x[[3]])), mf)
                 im <- as(ff, "sparseMatrix") # transpose of indicators
		 ## Could well be that we should rather check earlier .. :
		 if(!isTRUE(validObject(im, test=TRUE)))
		     stop("invalid conditioning factor in random effect: ", format(x[[3]]))

                 mm <- model.matrix(eval(substitute(~ expr, # model matrix
                                                    list(expr = x[[2]]))),
                                    mf)
                 if (rmInt) {
                     if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
                     if (ncol(mm) < 2)
                         stop("lhs of a random-effects term cannot be an intercept only")
                     mm <- mm[ , -icol , drop = FALSE]
                 }
                 list(f = ff,
                      Zt = drop0(do.call(rBind,
                      lapply(seq_len(ncol(mm)),
                             function(j) {im@x <- mm[,j]; im}))),
                      ST = matrix(0, ncol(mm), ncol(mm),
                      dimnames = list(colnames(mm), colnames(mm))))
             })
    dd <-
        VecFromNames(dimsNames, "integer",
                     c(list(n = nrow(mf), p = ncol(fr$X), nt = length(fl),
                            q = sum(sapply(fl, function(el) nrow(el$Zt)))),
                       dimsDefault))
    ## order terms by decreasing number of levels in the factor but don't
    ## change the order if this is already true
    nlev <- sapply(fl, function(el) length(levels(el$f)))
    ## determine the number of random effects at this point
    if (any(diff(nlev)) > 0) fl <- fl[rev(order(nlev))]
    ## separate the terms from the factor list
    trms <- lapply(fl, "[", -1)
    names(trms) <- NULL
    fl <- lapply(fl, "[[", "f")
    attr(fl, "assign") <- seq_along(fl)
    ## check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        ## check that the lengths of the number of levels coincide
        fl <- fl[match(ufn, fnms)]
        attr(fl, "assign") <- match(fnms, ufn)
    }
    names(fl) <- ufn
    ## check for nesting of factors
    dd["nest"] <- all(sapply(seq_along(fl)[-1],
                             function(i) isNested(fl[[i-1]], fl[[i]])))

    list(trms = trms, fl = fl, dims = dd)
}

mkZt <- function(fr, FL, start, s = 1L)
### Create the standard versions of flist, Zt, Gp, ST, A and L.
### Update dd.
{
    dd <- FL$dims
    fl <- FL$fl
    asgn <- attr(fl, "assign")
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn
    trms <- FL$trms
    Ztl <- lapply(trms, `[[`, "Zt")
    Zt <- do.call(rBind, Ztl)
    Zt@Dimnames <- vector("list", 2)
    ST <- new("ST", ST = lapply(trms, `[[`, "ST"),
              Gp = unname(c(0L, cumsum(sapply(Ztl, nrow)))))
    rm(Ztl, FL)                         # because they could be large
    .Call("ST_initialize", ST, Zt, PACKAGE = "lme4")
    A <- .Call("ST_create_A", ST, Zt, PACKAGE = "lme4")
    L <- Cholesky(tcrossprod(A), perm = TRUE, LDL = FALSE, Imult = 1)
    perm <- L@perm
    L@perm <- seq_len(nrow(A)) - 1L
    if (!is.null(start) && checkSTform(ST, start)) ST <- start
    
### FIXME: Check number of variance components versus number of
### levels in the factor for each term. Warn or stop as appropriate
    bds <- .Call("ST_bounds", ST, PACKAGE = "lme4")
    if (length(bds) %% 2) stop("length of bounds vector must be even")
    dd["np"] <- as.integer(length(bds)/2) # number of parameters in optimization
    X <- fr$X
    q <- nrow(Zt)
    y <- unname(fr$Y)
    n <- length(y)
    p <- ncol(X)
    wts <- unname(fr$wts)
    ans <- new(Class = "mer",
               env = new.env(),
               nlmodel = (~I(x))[[2]],
               frame = fr$mf,
               flist = fl,
               X = X,
               Zt = Zt,
               pWt = wts,
               offset = unname(fr$off),
               y = y,
               dims = dd,
               A = A,
               L = L,
               deviance = VecFromNames(devNames, "numeric"),
               fixef = fr$fixef,
               ranef = numeric(q),
               u = numeric(q),
               perm = perm,
               eta = numeric(n),
               mu = numeric(n),
               resid = numeric(n),
               sqrtrWt = sqrt(wts),
               RZX = matrix(0, q, p),
               RX = matrix(0, p, p),
               ghx = numeric(0),
               ghw = numeric(0))
    .Call("mer_update_mu", ans, PACKAGE = "lme4")
    list(mer = ans, ST = ST)
}

mer.tst <-
    function(formula, data, family = NULL, REML = TRUE,
             control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    mc <- match.call()
    if (!is.null(family)) {             # call glmer
        mc[[1]] <- as.name("glmer")
        return(eval.parent(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)
    fr <- lme4:::lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    FL <- lmerFactorList(formula, fr, 0L, 0L) # flist, Zt, dims
    Ztl <- mkZt(fr, FL, NULL)
    Ztl
}

fm1 <- mer.tst(Yield ~ (1|Batch), Dyestuff)
