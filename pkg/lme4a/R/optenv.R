#' Optimization-related methods for the environment class

if (FALSE) {
setMethod("getPars", representation(x = "environment"),
          function(x, ...)
      {
          if (is.null(active <- x$.active) || !is.character(active)
              || !all(active %in% objects(x)))
              stop("environment passed to getPars must contain a character object named .active")
          plist <- lapply(active, function(nm) getPars(x[[nm]], x, ...))
          if (is.null(x$.assign))
              x$.assign <- rep.int(seq_along(plist),
                                   unlist(lapply(plist, length)))
          unlist(plist)
      })

setMethod("getBounds", representation(x = "environment"),
          function(x, ...)
      {
          if (is.null(active <- x$.active) || !is.character(active)
              || !all(active %in% objects(x)))
              stop("environment passed to getBounds must contain a character object named .active")
          do.call(rbind, lapply(active, function(nm) getBounds(x[[nm]], x, ...)))
      })

setMethod("setPars", representation(x = "environment", pars = "numeric"),
          function(x, pars, ...)
      {
	  if (is.null(active <- x$.active) || !is.character(active)
	      || !all(active %in% objects(x)))
	      stop(gettextf("environment passed to setPars must contain a character object named .active"))
	  acseq <- seq_along(active)
	  if (is.null(assign <- x$.assign) || !is.integer(assign) ||
	      !all(assign %in% acseq))
	      stop(gettextf("object named .assign missing in environment passed to setPars or of incorrect form"))
	  if (length(assign) != length(pars))
	      stop(gettextf("'pars' vector is not of length %d",
			    length(assign)))

          ff <- function(i)
          {
              ans <- try(setPars(x[[ active[i] ]],pars[assign == i], x, ...),
                         silent = TRUE)
              if (class(ans) == "try-error") return(NA)
              ans
          }
          sum(unlist(lapply(acseq, ff)))
      })
}
setMethod("getPars", "optenv", function(x, ...) x@getPars())

setMethod("getBounds", "optenv", function(x, ...) x@getBounds())

setMethod("setPars", representation(x = "optenv", pars = "numeric"), function(x, pars, ...) x@setPars(pars))

setMethod("env", "optenv", function(x, ...) environment(x@getPars))
if (FALSE) {
setMethod("ranef", signature(object = "environment"),
          function(object, postVar = FALSE, drop = FALSE, whichel = names(ans), ...)
      {
          ## evaluate the list of matrices
          ml <- ranef(object$rCF, u = object$u, perm = object$perm,
                       postVar = postVar, ...)
          ## produce a list of data frames corresponding to factors, not terms
          fl <- object$flist
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
              pV <- .Call(ST_postVar, object$rCF, object$L,
                          object$perm, object$flist, whchL)
### need to multiply by sigma^2
              dd <- object$dims
              sc <- 1
              if (dd["useSc"])
                  sc <- object$deviance[if (dd["REML"]) "sigmaREML" else "sigmaML"]
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
}

