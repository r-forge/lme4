### Utilities for parsing the mixed model formula

##' From the right hand side of a formula for a mixed-effects model,
##' determine the pairs of expressions that are separated by the
##' vertical bar operator.
##'
##' @title Determine random-effects expressions from a formula
##' @seealso \code{\link{formula}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
##' @param term a mixed-model formula
##' @return pairs of expressions that were separated by vertical bars
##' @section Note: This function is called recursively on individual
##' terms in the model, which is why the argument is called \code{term} and not
##' a name like \code{form}, indicating a formula.
##' @examples
##' findbars(f1 <- Reaction ~ Days + (Days|Subject))
##' ## => list( Days | Subject )
##' findbars(y ~ Days + (1|Subject) + (0+Days|Subject))
##' ## => list of length 2:  list ( 1 | Subject ,  0+Days|Subject)
##' \dontshow{
##' stopifnot(identical(findbars(f1),
##'                     list(expression(Days | Subject)[[1]])))
##' }
##' @family utilities
##' @keywords models utilities
##' @export
findbars <- function(term)
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
    expandSlash(fb(term))
}

##' Remove the random-effects terms from a mixed-effects formula,
##' thereby producing the fixed-effects formula.
##'
##' @title Omit terms separated by vertical bars in a formula
##' @param term the right-hand side of a mixed-model formula
##' @return the fixed-effects part of the formula
##' @section Note: This function is called recursively on individual
##' terms in the model, which is why the argument is called \code{term} and not
##' a name like \code{form}, indicating a formula.
##' @examples
##' nobars(Reaction ~ Days + (Days|Subject)) ## => Reaction ~ Days
##' @seealso \code{\link{formula}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
##' @family utilities
##' @keywords models utilities
##' @export
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

##' Substitute the '+' function for the '|' function in a mixed-model
##' formula.  This provides a formula suitable for the current
##' model.frame function.
##'
##' @title "Sub[stitute] Bars"
##' @param term a mixed-model formula
##' @return the formula with all | operators replaced by +
##' @section Note: This function is called recursively on individual
##' terms in the model, which is why the argument is called \code{term} and not
##' a name like \code{form}, indicating a formula.
##' @examples
##' subbars(Reaction ~ Days + (Days|Subject)) ## => Reaction ~ Days + (Days + Subject)
##' @seealso \code{\link{formula}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
##' @family utilities
##' @keywords models utilities
##' @export
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

##' Does every level of f1 occur in conjunction with exactly one level
##' of f2? The function is based on converting a triplet sparse matrix
##' to a compressed column-oriented form in which the nesting can be
##' quickly evaluated.
##'
##' @title Is f1 nested within f2?
##'
##' @param f1 factor 1
##' @param f2 factor 2
##'
##' @return TRUE if factor 1 is nested within factor 2
##' @examples
##' with(Pastes, isNested(cask, batch))   ## => FALSE
##' with(Pastes, isNested(sample, batch))  ## => TRUE
##' @export
isNested <- function(f1, f2)
{
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    stopifnot(length(f1) == length(f2))
    k <- length(levels(f1))
    sm <- as(new("ngTMatrix",
		 i = as.integer(f2) - 1L,
		 j = as.integer(f1) - 1L,
		 Dim = c(length(levels(f2)), k)),
             "CsparseMatrix")
    all(sm@p[2:(k+1L)] - sm@p[1:k] <= 1L)
}

subnms <- function(form, nms) {
    ## Recursive function applied to individual terms
    sbnm <- function(term)
    {
        if (is.name(term))
            if (any(term == nms)) return(0) else return(term)
        switch(length(term),
               return(term),
           {
               term[[2]] <- sbnm(term[[2]])
               return(term)
           },
           {
               term[[2]] <- sbnm(term[[2]])
               term[[3]] <- sbnm(term[[3]])
               return(term)
           })
        NULL
    }
    sbnm(form)
}

## Check for a constant term (a literal 1) in an expression
##
## In the mixed-effects part of a nonlinear model formula, a constant
## term is not meaningful because every term must be relative to a
## nonlinear model parameter.  This function recursively checks the
## expressions in the formula for a a constant, calling stop() if
## such a term is encountered.
## @title Check for constant terms.
## @param expr an expression
## @return NULL.  The function is executed for its side effect.
chck1 <- function(expr) {
    if ((le <- length(expr)) == 1) {
        if (is.numeric(expr) && expr == 1) 
            stop("1 is not meaningful in a nonlinear model formula")
        return()
    } else
    for (j in seq_len(le)[-1]) Recall(expr[[j]])
}

##' Check and manipulate the formula for a nonlinear model.
##'
##' The model formula for a nonlinear mixed-effects model is of the form
##' \code{resp ~ nlmod ~ mixed} where \code{"resp"} is an expression
##' (usually just a name) for the response, \code{nlmod} is the
##' call to the nonlinear model function, and \code{mixed} is the
##' mixed-effects formula defining the linear predictor for the
##' parameter matrix.  If the formula is to be used for optimizing
##' designs, the \code{"resp"} part can be omitted.
##'
##' 
##' @title Manipulate a nonlinear model formula.
##' @param form The model formula
##' @param pnames a character vector of parameter names.
##' @param need3 logical scalar indicating if the model formula must
##'      be a three-part formula.
##' @return a list with components
##'  \item{"nlmod"}{the call to the nonlinear model function}
##'  \item{"fr.form"}{the formula for the call to \code{\link{model.frame}}}
##'  \item{"fe"}{the fixed-effects model formula}
##'  \item{"re"}{the list of random-effects terms}
##' @export
##' @family utilities
nlformula <- function(form, pnames, need3 = TRUE) {
    form <- as.formula(form)
    if (length(form) < 3L)
        stop("formula must include nonlinear and mixed-effects terms")
    chck1(meform <- form[[3L]])
    nlform <- as.formula(form[[2L]])
    environment(nlform) <- environment(form)
    if ((lnl <- length(nlform)) < 3L && need3) stop("formula must be a 3-part formula")

    pnameexpr <- parse(text=paste(pnames, collapse='+'))[[1]]
    if (length(nb <- nobars(meform)))
        nb <- substitute(~ 0 + nb + pnameexpr)
    else 
        nb <- substitute(~ 0 + pnameexpr)
    fb <- lapply(findbars(meform),
                 function(expr) {expr[[2]] = substitute(0+foo, list(foo=expr[[2]]));expr})
    frexpr <- parse(text= paste(setdiff(all.vars(form), pnames), collapse=' + '))[[1]]
    fr.form <- if (lnl < 3L) substitute(~ frexpr) else
    substitute(resp ~ frexpr, list(resp = nlform[[2]], frexpr=frexpr))
    environment(fr.form) <- environment(form)
    list(nlmod = as.call(nlform[[length(nlform)]]), fr.form=fr.form, fe=nb, re=fb)
}
