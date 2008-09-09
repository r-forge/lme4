## Preliminary version of a function to fit a model with carryover for
## the effect of one grouping factor (e.g. teacher) on another
## grouping factor (e.g. student)

lagged_factor <- function(fr, idvar = "id", timevar = "Time", factors = "inst", sep = ".")
{
    ## drop unused levels of any factors in fr
    fr <- do.call(data.frame, lapply(fr, "[", drop = TRUE))
    ## extract names and check arguments for consistency
    fnms <- names(fr)
    idvar <- as.character(idvar)
    timevar <- as.character(timevar)
    factors <- as.character(factors)
    stopifnot(length(idvar) == 1, length(timevar) == 1,
              idvar %in% fnms, timevar %in% fnms,
              all(factors %in% fnms),
              all(unlist(lapply(fr[factors], is.factor))))
    v.names <- setdiff(fnms, c(idvar, timevar))
    wide <- reshape(fr, direction = "wide", idvar = idvar, timevar = timevar, sep = sep,
                    v.names = v.names)
    nts <- seq_len(nt <- length(tlevs <- levels(fr[[timevar]])))
    varying <- attr(wide, "reshapeWide")$varying
    for (nm in factors) {
        levs <- levels(fr[[nm]])
        vnms <- paste(nm, tlevs, sep = sep)
        ff <- cbind(wide[, vnms],
                    data.frame(NAs = factor(rep(NA, nrow(wide)),
                               levels = levs))[, rep.int(1L, nt - 1)])
        arr <- array(unlist(ff), c(nrow(wide), nt + 1L, nt))[, nts, nts]
        for (face in nts[-1]) {
            nms <- paste(nm, face, sep, tlevs, sep = '')
            for (j in nts) wide[[nms[j]]] <- factor(arr[,j,face], levels = levs)
            varying <- rbind(varying, nms)
        }
    }
    attr(wide, "reshapeWide")$varying <- varying
    attr(wide, "reshapeWide")$v.names <- c(v.names, as.vector(outer(nm, nts[-1], paste, sep = sep)))
    reshape(wide)
}

carryOver <-
    function(formula, data, family = NULL, REML = TRUE, nAGQ = 1,
             control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, idvar = "id", timevar = "Time",
             factors = "sch", ...)
{
    lmerc <- mc <- match.call()
    ## Call lmer without Znew and pre and with doFit = FALSE
    lmerc$idvar <- lmerc$timevar <- lmerc$factors <- NULL
    mc$doFit <- FALSE
    mc[[1]] <- as.name("lmer")
    lf <- eval.parent(mc)
browser()
    ans <- do.call(if (!is.null(lf$glmFit))
                   lme4:::glmer_finalize else lme4:::lmer_finalize, lf)
    ans@call <- match.call()
    ans
}
