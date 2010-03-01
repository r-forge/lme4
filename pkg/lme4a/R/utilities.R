### Utilities for parsing the mixed model formula

##' Determine random-effects expressions from a formula

##' From the right hand side of a formula for a mixed-effects model,
##' determine the pairs of expressions that are separated by the
##' vertical bar operator.

##' @param term the right-hand side of a mixed-model formula

##' @return pairs of expressions that were separated by vertical bars
findbars <- function(expr)
{
    ## Recursive function applied to individual terms
    fb <- function(term)
    {
        if (is.name(term) || !is.language(term)) return(NULL)
        if (term[[1]] == as.name("(")) return(fb(term[[2]]))
        stopifnot(is.call(term))
        if (term[[1]] == as.name('|')) return(term)
        if (length(term) == 2) return(fb(term[[2]]))
        c(fb(term[[2]]), fb(term[[3]]))
    }
    ## Expand any slashes in the grouping factors returned by fb
    expandSlash <- function(bb)
    {
        ## Create the interaction terms for nested effects
        makeInteraction <- function(x)
        {
            if (length(x) < 2) return(x)
            trm1 <- makeInteraction(x[[1]])
            trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
            list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
        }
        ## Return the list of '/'-separated terms
        slashTerms <- function(x)
        {
            if (!("/" %in% all.names(x))) return(x)
            if (x[[1]] != as.name("/"))
                stop("unparseable formula for grouping factor")
            list(slashTerms(x[[2]]), slashTerms(x[[3]]))
        }
        
        if (!is.list(bb)) return(expandSlash(list(bb)))
        ## lapply(unlist(... - unlist returns a flattened list
        unlist(lapply(bb, function(x) {
            if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
                return(lapply(unlist(makeInteraction(trms)),
                              function(trm) substitute(foo|bar,
                                                       list(foo = x[[2]],
                                                            bar = trm))))
            x
        }))
    }
    expandSlash(fb(expr))
}

##' Omit terms separated by vertical bars in a formula
##'
##' Remove the random-effects terms from a mixed-effects formula,
##' thereby producing the fixed-effects formula.
##' @param term the right-hand side of a mixed-model formula

##' @return the fixed-effects part of the formula
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

##' "Sub[stitute] Bars"

##' Substitute the '+' function for the '|' function in a mixed-model
##' formula.  This provides a formula suitable for the current
##' model.frame function.
##' @param term the right-hand side of a mixed-model formula

##' @return the formula with all | operators replaced by +
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

##' Expand slashes in grouping factors

##' Expand any slashes in the grouping factors returned by fb
##' @param bb a list of terms returned by fb

##' @return a list of pairs of expressions that generate
##'    random-effects terms with slashes expanded.
expandSlash <- function(bb)
{
    ## Create the interaction terms for nested effects
    makeInteraction <- function(x)
    {
        if (length(x) < 2) return(x)
        trm1 <- makeInteraction(x[[1]])
        trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
        list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
    }
    ## Return the list of '/'-separated terms
    slashTerms <- function(x)
    {
        if (!("/" %in% all.names(x))) return(x)
        if (x[[1]] != as.name("/"))
            stop("unparseable formula for grouping factor")
        list(slashTerms(x[[2]]), slashTerms(x[[3]]))
    }

    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## lapply(unlist(... - unlist returns a flattened list
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
