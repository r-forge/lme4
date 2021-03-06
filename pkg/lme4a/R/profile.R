setGeneric("dropX",
           function(x, which, fw) standardGeneric("dropX"))
## This is a hack.  The preferred approach is to write a
## subset method for the ddenseModelMatrix and dsparseModelMatrix
## classes

.modelMatrixDrop <- function(mm, w) {
    ll <- list(Class = class(mm),
               assign = mm@assign[-w],
               contrasts = mm@contrasts)
    X <- mm[, -w, drop = FALSE]
    ll <- c(ll, lapply(structure(slotNames(X), .Names=slotNames(X)),
                       function(nm) slot(X, nm)))
    do.call("new", ll)
}

##' Drop the which'th column from the fixed-effects model matrix in
##' fe, the deFeMod slot.  Add fw times the dropped column to
##' resp@offset and store values to implement
##' profiling of the fixed-effects parameters.
##'
##' @title Drop the which'th column from X in an merMod object
##'
##' @param obj
##' @param which the column to drop.  Must have 1 <= which <= ncol(obj@fe@X)
##' @param fw the value of the which'th fixed effect which will
##'     be held constant.
##' @return a revised merMod object
setMethod("dropX", "merMod",
          function(x, which, fw)
      {
          w <- as.integer(which)[1]
          fw <- as.numeric(fw)[1]
          fe <- x@fe
          resp <- x@resp
          p <- length(fe@coef)
          stopifnot(0 < w, w <= p)
          Xw <-fe@X[, w, drop = TRUE]
          X <- .modelMatrixDrop(fe@X, w)
          VtV <- crossprod(X)
          ll <- list(Class = class(fe),
                     coef = fe@coef[-w],
                     RZX = fe@RZX[, -w, drop = FALSE],
                     fac = chol(VtV),
                     X = X)
          ## offset calculated from fixed parameter value
          resp@offset <- Xw * fw + resp@offset
          newfe <- do.call("new", ll)
          new("merMod", re = x@re, fe = newfe, resp = resp,
              call = x@call, devcomp = x@devcomp, frame = x@frame)
      })


##' extract only the y component from a prediction
predy <- function(sp, vv) predict(sp, vv)$y

stripExpr <- function(ll, nms) {
    stopifnot(inherits(ll, "list"), is.character(nms))
    lsigNm <- which(nms == ".lsig")
    sigNms <- grep("^.sig", nms)
    sigsub <- as.integer(substring(nms[sigNms], 5))
    fLevs <- as.expression(nms)
    fLevs[lsigNm] <- expression(log(sigma))
    fLevs[sigNms] <- parse(text=paste("sigma[", sigsub, "]"))
    levsExpr <- substitute(strip.custom(factor.levels=foo), list(foo=fLevs))
    llNms <- names(ll)
    snames <- c("strip", "strip.left")
    if (all(!(snames %in% llNms))) {
        ll$strip <- levsExpr
    } else {
        lapply(snames, function(nm) {
            if (nm %in% llNms) {
                vv <- ll[[nm]]
                if (is.logical(vv) && vv) ll[[nm]] <<- levsExpr
            }
        })
    }
    ll
}

## A lattice-based plot method for profile objects
xyplot.thpr <-
    function (x, data = NULL,
              levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), 1)),
              conf = c(50, 80, 90, 95, 99)/100,
              absVal = FALSE, ...)
{
    levels <- sort(levels[is.finite(levels) & levels > 0])
    spl <- attr(x, "forward")
    bspl <- attr(x, "backward")
    zeta <- c(-rev(levels), 0, levels)
    fr <- data.frame(zeta = rep.int(zeta, length(spl)),
                     pval = unlist(lapply(bspl, predy, zeta)),
                     pnm = gl(length(spl), length(zeta), labels = names(spl)))
    if (length(ind <- which(is.na(fr$pval)))) {
        fr[ind, "zeta"] <- 0
        for (i in ind)
### FIXME: Should check which bound has been violated, although it
### will almost always be the minimum.
            fr[i, "pval"] <- min(spl[[fr[i, "pnm"] ]]$knots)
    }
    ylab <- expression(zeta)
    if (absVal) {
        fr$zeta <- abs(fr$zeta)
        ylab <- expression("|" * zeta * "|")
    }
    ll <- c(list(...),
            list(x = zeta ~ pval | pnm, data=fr,
                 scales = list(x = list(relation = 'free')),
                 ylab = ylab, xlab = NULL,
                 panel = function(x, y, ...)
             {
                 panel.grid(h = -1, v = -1)
                 lsegments(x, y, x, 0, ...)
                 lims <- current.panel.limits()$xlim
                 myspl <- spl[[panel.number()]]
                 krange <- range(myspl$knots)
                 pr <- predict(myspl,
                               seq(max(lims[1], krange[1]),
                                   min(lims[2], krange[2]), len = 101))
                 if (absVal) {
                     pr$y <- abs(pr$y)
                     y[y == 0] <- NA
                     lsegments(x, y, rev(x), y)
                 } else {
                     panel.abline(h = 0, ...)
                 }
                 panel.lines(pr$x, pr$y)
             }))
    do.call(xyplot, stripExpr(ll, names(spl)))
}

confint.thpr <- function(object, parm, level = 0.95, zeta, ...)
{
    bak <- attr(object, "backward")
    bnms <- names(bak)
    if (missing(parm)) parm <- bnms
    else if (is.numeric(parm)) parm <- bnms[parm]
    parm <- intersect(as.character(parm), bnms)
    cn <- NULL
    if (missing(zeta)) {
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        zeta <- qnorm(a)
        cn <- stats:::format.perc(a, 3)
    }
    ci <- t(sapply(parm, function(nm) predy(bak[[nm]], zeta)))
    colnames(ci) <- cn
    ci
}

##' Convert x-cosine and y-cosine to average and difference.

##' Convert the x-cosine and the y-cosine to an average and difference
##' ensuring that the difference is positive by flipping signs if
##' necessary
ad <- function(xc, yc)
{
    a <- (xc + yc)/2
    d <- (xc - yc)
    cbind(ifelse(d > 0, a, -a), abs(d))
}

##' convert d versus a (as an xyVector) and level to a matrix of taui and tauj
tauij <- function(xy, lev) lev * cos(xy$x + outer(xy$y/2, c(-1, 1)))

## safe arc-cosine
sacos <- function(x) acos(pmax.int(-0.999, pmin.int(0.999, x)))

cont <- function(sij, sji, levels, nseg = 101)
{
    ada <- array(0, c(length(levels), 2, 4))
    ada[, , 1] <- ad(0, sacos(predy(sij,  levels)/levels))
    ada[, , 2] <- ad(sacos(predy(sji, levels)/levels), 0)
    ada[, , 3] <- ad(pi, sacos(predy(sij, -levels)/levels))
    ada[, , 4] <- ad(sacos(predy(sji, -levels)/levels), pi)
    pts <- array(0, c(length(levels), nseg + 1, 2))
    for (i in seq_along(levels))
        pts[i, ,] <- tauij(predict(periodicSpline(ada[i, 1, ], ada[i, 2, ]),
                                   nseg = nseg), levels[i])
    levs <- c(-rev(levels), 0, levels)
    list(tki = predict(sij, levs), tkj = predict(sji, levs), pts = pts)
}

## Contours are for the marginal two-dimensional regions (i.e. using
## df = 2)
splom.thpr <-
    function (x, data, ## unused - only for compatibility with generic
              levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), 2)),
              conf = c(50, 80, 90, 95, 99)/100, ...)
{
    mlev <- max(levels)
    spl <- attr(x, "forward")
    frange <- sapply(spl, function(x) range(x$knots))
    bsp <- attr(x, "backward")
    brange <- sapply(bsp, function(x) range(x$knots))
    pfr <- do.call(cbind, lapply(bsp, predy, c(-mlev, mlev)))
    pfr[1, ] <- pmax.int(pfr[1, ], frange[1, ], na.rm = TRUE)
    pfr[2, ] <- pmin.int(pfr[2, ], frange[2, ], na.rm = TRUE)
    nms <- names(spl)
    ## Create data frame fr of par. vals in zeta coordinates
    fr <- x[, -1]
    for (nm in nms) fr[[nm]] <- predy(spl[[nm]], fr[[nm]])
    fr1 <- fr[1, nms]
    ## create a list of lists with the names of the parameters
    traces <- lapply(fr1, function(el) lapply(fr1, function(el1) list()))
    for (j in seq_along(nms)[-1]) {
        for (i in seq_len(j - 1)) {
            fri <- subset(fr, .par == nms[i])
            sij <- interpSpline(fri[ , i], fri[ , j])
            frj <- subset(fr, .par == nms[j])
            sji <- interpSpline(frj[ , j], frj[ , i])
            ll <- cont(sij, sji, levels)
            traces[[j]][[i]] <- list(sij = sij, sji = sji, ll = ll)
        }
    }
    ## panel function for lower triangle
    lp <- function(x, y, groups, subscripts, i, j, ...) {
        tr <- traces[[j]][[i]]
        grid::pushViewport(viewport(xscale = c(-1.07, 1.07) * mlev,
                              yscale = c(-1.07, 1.07) * mlev))
        dd <- sapply(current.panel.limits(), diff)/50
        psij <- predict(tr$sij)
        ll <- tr$ll
        ## now do the actual plotting
        panel.grid(h = -1, v = -1)
        llines(psij$y, psij$x, ...)
        llines(predict(tr$sji), ...)
        with(ll$tki, lsegments(y - dd[1], x, y + dd[1], x, ...))
        with(ll$tkj, lsegments(x, y - dd[2], x, y + dd[2], ...))
        for (k in seq_along(levels)) llines(ll$pts[k, , ], ...)
        grid::popViewport(1)
    }
    ## panel function for upper triangle
    up <- function(x, y, groups, subscripts, i, j, ...) {
        ## panels are transposed so reverse i and j
        jj <- i
        ii <- j
        tr <- traces[[jj]][[ii]]
        ll <- tr$ll
        pts <- ll$pts
        limits <- current.panel.limits()
        psij <- predict(tr$sij)
        psji <- predict(tr$sji)
        ## do the actual plotting
        panel.grid(h = -1, v = -1)
        llines(predy(bsp[[ii]], psij$y), predy(bsp[[jj]], psij$x), ...)
        llines(predy(bsp[[ii]], psji$x), predy(bsp[[jj]], psji$y), ...)
        for (k in seq_along(levels))
            llines(predy(bsp[[ii]], pts[k, , 1]),
                   predy(bsp[[jj]], pts[k, , 2]), ...)
    }
    dp <- function(x = NULL,            # diagonal panel
                   varname = NULL, limits, at = NULL, lab = NULL,
                   draw = TRUE,

                   varname.col = add.text$col,
                   varname.cex = add.text$cex,
                   varname.lineheight = add.text$lineheight,
                   varname.font = add.text$font,
                   varname.fontfamily = add.text$fontfamily,
                   varname.fontface = add.text$fontface,

                   axis.text.col = axis.text$col,
                   axis.text.alpha = axis.text$alpha,
                   axis.text.cex = axis.text$cex,
                   axis.text.font = axis.text$font,
                   axis.text.fontfamily = axis.text$fontfamily,
                   axis.text.fontface = axis.text$fontface,

                   axis.line.col = axis.line$col,
                   axis.line.alpha = axis.line$alpha,
                   axis.line.lty = axis.line$lty,
                   axis.line.lwd = axis.line$lwd,
                   i, j, 
                   ...)
    {
        n.var <- eval.parent(expression(n.var))
        add.text <- trellis.par.get("add.text")
        axis.line <- trellis.par.get("axis.line")
        axis.text <- trellis.par.get("axis.text")
        if (!is.null(varname))
            grid::grid.text(varname,
                            gp =
                            gpar(col = varname.col,
                                 cex = varname.cex,
                                 lineheight = varname.lineheight,
                                 fontface = lattice:::chooseFace(varname.fontface,
                                 varname.font),
                                 fontfamily = varname.fontfamily))
        if (draw) {
            at <- pretty(limits)
            sides <- c("left", "top")
            if (j == 1) sides <- "top"
            if (j == n.var) sides <- "left"
            for (side in sides)
                panel.axis(side = side,
                           at = at,
                           labels = format(at, trim = TRUE),
                           tick = TRUE,
                           check.overlap = TRUE,
                           half = side == "top" && j > 1,

                           tck = 1, rot = 0,

                           text.col = axis.text.col,
                           text.alpha = axis.text.alpha,
                           text.cex = axis.text.cex,
                           text.font = axis.text.font,
                           text.fontfamily = axis.text.fontfamily,
                           text.fontface = axis.text.fontface,

                           line.col = axis.line.col,
                           line.alpha = axis.line.alpha,
                           line.lty = axis.line.lty,
                           line.lwd = axis.line.lwd)
            lims <- c(-1.07, 1.07) * mlev
            grid::pushViewport(viewport(xscale = lims, yscale = lims))
            side <- ifelse(j == 1, "right", "bottom")
            which.half <- ifelse(j == 1, "lower", "upper")
            at <- pretty(lims)
            panel.axis(side = side, at = at, labels = format(at, trim = TRUE),
                       tick = TRUE, half = TRUE, which.half = which.half,
                       tck = 1, rot = 0,

                       text.col = axis.text.col,
                       text.alpha = axis.text.alpha,
                       text.cex = axis.text.cex,
                       text.font = axis.text.font,
                       text.fontfamily = axis.text.fontfamily,
                       text.fontface = axis.text.fontface,

                       line.col = axis.line.col,
                       line.alpha = axis.line.alpha,
                       line.lty = axis.line.lty,
                       line.lwd = axis.line.lwd)
            grid::popViewport(1)
        }
    }

    splom(~ pfr, lower.panel = lp, upper.panel = up, diag.panel = dp, ...)
}

##' Transform an lmer profile to the scale of the logarithm of the
##' standard deviation of the random effects.
##' @title Transform an lmer profile to the logarithm scale
##' @param x an object that inherits from class "thpr"
##' @param base the base of the logarithm.  Defaults to natural
##'        logarithms
##'
##' @return an lmer profile like x with all the .sigNN parameters
##'      replaced by .lsigNN.  The forward and backward splines for
##'      these parameters are recalculated.
log.thpr <- function (x, base = exp(1))
{
    cn <- colnames(x)
    sigs <- grep("^\\.sig[0-9][0-9]", cn)
    if (length(sigs)) {
        colnames(x) <- sub("^\\.sig", ".lsig", cn)
        levels(x$.par) <- sub("^\\.sig", ".lsig", levels(x$.par))
        names(attr(x, "backward")) <-
            names(attr(x, "forward")) <-
                sub("^\\.sig", ".lsig", names(attr(x, "forward")))
        for (nm in colnames(x)[sigs]) {
            x[[nm]] <- log(x[[nm]], base = base)
            fr <- subset(x, .par == nm & is.finite(x[[nm]]))
            form <- eval(substitute(.zeta ~ nm, list(nm = as.name(nm))))
            attr(x, "forward")[[nm]] <- interpSpline(form, fr)
            attr(x, "backward")[[nm]] <- backSpline(attr(x, "forward")[[nm]])
        }
        ## eliminate rows the produced non-finite logs
        x <- x[apply(is.finite(as.matrix(x[, sigs])), 1, all),]
    }
    x
}

## FIXME - change this so it uses .sigsq too
varpr <- function (x)
{
    cn <- colnames(x)
    sigs <- grep("^\\.sig[0-9][0-9]", cn)
    if (length(sigs)) {
        colnames(x) <- sub("^\\.sig", ".sigsq", cn)
        levels(x$.par) <- sub("^\\.sig", ".sigsq", levels(x$.par))
        names(attr(x, "backward")) <-
            names(attr(x, "forward")) <-
                sub("^\\.sig", ".sigsq", names(attr(x, "forward")))
        for (nm in colnames(x)[sigs]) {
            x[[nm]] <- x[[nm]]^2
            fr <- subset(x, .par == nm & is.finite(x[[nm]]))
            form <- eval(substitute(.zeta ~ nm, list(nm = as.name(nm))))
            attr(x, "forward")[[nm]] <- interpSpline(form, fr)
            attr(x, "backward")[[nm]] <- backSpline(attr(x, "forward")[[nm]])
        }
        ## eliminate rows the produced non-finite logs
        x <- x[apply(is.finite(as.matrix(x[, sigs])), 1, all),]
    }
    x
}

##' The deviance is profiled with respect to the fixed-effects
##' parameters but not with respect to sigma, which is expressed
##' on the logarithmic scale, lsigma. The other parameters are on
##' the standard deviation scale, not the theta scale.
##'
##' @title Return a function for evaluation of the deviance.
##' @param fm a fitted model of class merMod
##' @return a function for evaluating the deviance in the extended
##'     parameterization.  This is profiled with respect to the
##'     fixed-effects but not with respect to sigma.
devfun2 <- function(fm)
{
    stopifnot(is(fm, "merMod"), is(fm@resp, "lmerResp"))
    fm1 <- refitML(fm)
    th <- fm1@re@theta
    lth <- length(th)
    rm(fm)
    basedev <- unname(deviance(fm1))
    sig <- sigma(fm1)
    lsig <- log(sig)
    np <- lth + 1L
    ans <- function(pars)
    {
        stopifnot(is.numeric(pars), length(pars) == np)
        ## Assumption:  1) last parameter = log(sigma); 2) other pars are on SD-scale
        sigma <- exp(pars[np])
        sigsq <- sigma^2
        pp <- pars[-np]/sigma
        fv <- .Call(merDeviance, fm1, pp, fm1@fe@coef, fm1@re@u, 0L, 1L)
        aa <- attr(fv,"ldL2") + (attr(fv,"wrss")+attr(fv,"ussq"))/sigsq +
            length(fm1@resp@y) * log(2 * pi * sigsq)
        mostattributes(aa) <- attributes(fv)
        aa
    }
    opt <- c(sig * th, lsig)
    names(opt) <- c(sprintf(".sig%02d", seq_len(lth)), ".lsig")
    attr(ans, "optimum") <- c(opt, fixef(fm1)) # w/ names()
    attr(ans, "basedev") <- basedev
    attr(ans, "thopt") <- th
    attr(ans, "stderr") <- sig * sqrt(unscaledVar(RX = fm1@fe@fac))
    class(ans) <- "devfun"
    ans
}

profile.merMod <- function(fitted, alphamax = 0.01, maxpts = 100, delta = cutoff/8,
                           tr = 0, ...)
{
    dd <- devfun2(fitted)
    
    base <- attr(dd, "basedev")
    thopt <- attr(dd, "thopt")
    stderr <- attr(dd, "stderr")
    fm1 <- environment(dd)$fm1
    X.orig <- fm1@fe@X
    n <- length(fm1@resp@y)
    p <- length(fm1@fe@coef)
    
    ans <- lapply(opt <- attr(dd, "optimum"), function(el) NULL)
    bakspl <- forspl <- ans
    
    nptot <- length(opt)
    nvp <- nptot - p    # number of variance-covariance pars
    fe.orig <- opt[-seq_len(nvp)]
    res <- c(.zeta = 0, opt)
    res <- matrix(res, nr = maxpts, nc = length(res),
                  dimnames = list(NULL, names(res)), byrow = TRUE)
    cutoff <- sqrt(qchisq(1 - alphamax, nptot))
    
    ## helper functions
    
    ## nextpar calculates the next value of the parameter being
    ## profiled based on the desired step in the profile zeta
    ## (absstep) and the values of zeta and column cc for rows
    ## r-1 and r.  The parameter may not be below lower
    nextpar <- function(mat, cc, r, absstep, lower = -Inf) {
        rows <- r - (1:0)         # previous two row numbers
        pvals <- mat[rows, cc]
        zeta <- mat[rows, ".zeta"]
        if (!(denom <- diff(zeta)))
            stop("Last two rows have identical .zeta values")
        num <- diff(pvals)
        max(lower, pvals[2] + sign(num) * absstep * num / denom)
    }
    
    ## mkpar generates the parameter vector of theta and
    ## log(sigma) from the values being profiled in position w
    mkpar <- function(np, w, pw, pmw) {
        par <- numeric(np)
        par[w] <- pw
        par[-w] <- pmw
        par
    }
    
    ## fillmat fills the third and subsequent rows of the matrix
    ## using nextpar and zeta
### FIXME:  add code to evaluate more rows near the minimum if that
###        constraint was active.
    fillmat <- function(mat, lowcut, zetafun, cc) {
        nr <- nrow(mat)
        i <- 2L
        while (i < nr && abs(mat[i, ".zeta"]) <= cutoff &&
               mat[i, cc] > lowcut) {
            mat[i + 1L, ] <- zetafun(nextpar(mat, cc, i, delta, lowcut))
            i <- i + 1L
        }
        mat
    }
    
    lower <- c(fm1@re@lower, rep.int(-Inf, p + 1L ))
    seqnvp <- seq_len(nvp)
    lowvp <- lower[seqnvp]
    form <- .zeta ~ foo           # pattern for interpSpline formula
    
    for (w in seqnvp) {
        wp1 <- w + 1L
        start <- opt[seqnvp][-w]
        pw <- opt[w]
        lowcut <- lower[w]
        zeta <- function(xx) {
            ores <- bobyqa(start,
                           function(x) dd(mkpar(nvp, w, xx, x)),
                           lower = lowvp[-w])
            zz <- sign(xx - pw) * sqrt(ores$fval - base)
            c(zz, mkpar(nvp, w, xx, ores$par), attr(ores$fval, "beta"))
        }
        
### FIXME: The starting values for the conditional optimization should
### be determined from recent starting values, not always the global
### optimum values.
        
### Can do this by taking the change in the other parameter values at
### the two last points and extrapolating.
        
        ## intermediate storage for pos. and neg. increments
        pres <- nres <- res
        ## assign one row, determined by inc. sign, from a small shift
        nres[1, ] <- pres[2, ] <- zeta(pw * 1.01)
        ## fill in the rest of the arrays and collapse them
        bres <-
            as.data.frame(unique(rbind2(fillmat(pres,lowcut, zeta, wp1),
                                        fillmat(nres,lowcut, zeta, wp1))))
        bres$.par <- names(opt)[w]
        ans[[w]] <- bres[order(bres[, wp1]), ]
        form[[3]] <- as.name(names(opt)[w])
        
        bakspl[[w]] <- backSpline(forspl[[w]] <- interpSpline(form, bres))
    }
    
    offset <- fm1@resp@offset
    X <- fm1@fe@X
    
    for (j in seq_len(p)) {
        pres <-            # intermediate results for pos. incr.
            nres <- res    # and negative increments
        est <- opt[nvp + j]
        std <- stderr[j]
        Xw <- X[ , j, drop = TRUE]
        fmm1 <- dropX(fm1, j, est)
        fe.zeta <- function(fw) {
            fmm1@resp@offset <- Xw * fw + offset
            ores <- bobyqa(thopt, mkdevfun(fmm1),
                           lower = fmm1@re@lower)
            fv <- ores$fval
            sig <- sqrt((attr(fv, "wrss") + attr(fv, "ussq"))/length(Xw))
            c(sign(fw - est) * sqrt(ores$fval - base),
              ores$par * sig, log(sig), mkpar(p, j, fw, attr(fv, "beta")))
        }
        nres[1, ] <- pres[2, ] <- fe.zeta(est + delta * std)
        pp <- nvp + 1L + j
        bres <-
            as.data.frame(unique(rbind2(fillmat(pres,-Inf, fe.zeta, pp),
                                        fillmat(nres,-Inf, fe.zeta, pp))))
        thisnm <- names(fe.orig)[j]
        bres$.par <- thisnm
        ans[[thisnm]] <- bres[order(bres[, pp]), ]
        form[[3]] <- as.name(thisnm)
        bakspl[[thisnm]] <-
            backSpline(forspl[[thisnm]] <- interpSpline(form, bres))
    }
    
    ans <- do.call(rbind, ans)
    ans$.par <- factor(ans$.par)
    attr(ans, "forward") <- forspl
    attr(ans, "backward") <- bakspl
    row.names(ans) <- NULL
    class(ans) <- c("thpr", "data.frame")
    ans
}

dens <- function(pr, npts=201, upper=0.999) {
    npts <- as.integer(npts)
    stopifnot(inherits(pr, "thpr"), npts > 0,
              is.numeric(upper), 0.5 < upper, upper < 1)
    spl <- attr(pr, "forward")
    bspl <- attr(pr, "backward")
    zeta <- c(qnorm(1-upper), qnorm(upper))
    rng <- lapply(bspl, function(spl)
              {
                  rng <- predy(spl, zeta)
                  if (is.na(rng[1])) rng[1] <- 0
                  seq(rng[1], rng[2], len=npts)
              })
    fr <- data.frame(pval=unlist(rng),
                     pnm=gl(length(rng), npts, labels=names(rng)))
    dd <- list()
    for (nm in names(rng)) {
        zz <- predy(spl[[nm]], rng[[nm]])
        dd[[nm]] <- dnorm(zz) * predict(spl[[nm]], rng[[nm]], deriv=1)$y
    }
    fr$density <- unlist(dd)
    fr
}

densityplot.thpr <- function(x, data, ...) {
    ll <- c(list(...),
            list(x=density ~ pval|pnm,
                 data=dens(x),
                 type=c("l","g"),
                 scales=list(relation="free"),
                 xlab=NULL))
    do.call(xyplot, stripExpr(ll, names(attr(x, "forward"))))
}

varianceProf <- function(pr) {
    stopifnot(inherits(pr, "thpr"))
    spl <- attr(pr, "forward")
    onms <- names(spl)                  # names of original variables
    vc <- onms[grep("^.[l]*sig", onms)]  # variance components
    ans <- subset(pr, .par %in% vc, select=c(".zeta", vc, ".par"))
    ans$.par <- factor(ans$.par)        # drop unused levels
    if (".lsig" %in% vc) ans$.lsig <- exp(ans$.lsig)
    attr(ans, "forward") <- attr(ans, "backward") <- list()
    for (nm in vc) {
        ans[[nm]] <- ans[[nm]]^2
        fr <- subset(ans, .par == nm & is.finite(ans[[nm]]))
        form <- eval(substitute(.zeta ~ nm, list(nm = as.name(nm))))
        attr(ans, "forward")[[nm]] <- interpSpline(form, fr)
        attr(ans, "backward")[[nm]] <- backSpline(attr(ans, "forward")[[nm]])
    }
    ans
}
