#' Optimization-related methods for the environment class

setMethod("getPars", representation(x = "environment"),
          function(x, ...)
      {
          if (is.null(active <- x$.active) || !is.character(active)
              || !all(active %in% objects(x)))
              stop("environment passed to getPars must contain a character object named .active")
          plist <- lapply(active, function(nm) getPars(x[[nm]], x, ...))
          if (is.null(x$.assign))
              x$.assign <- rep.int(seq_along(plist), unlist(lapply(plist, length)))
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
              !all(assign %in% acseq) || length(assign) != length(pars))
              stop(gettextf("object named .assign missing in environment passed to setPars or of incorrect form"))

          sum(unlist(lapply(acseq, function(i) setPars(x[[ active[i] ]],
                                                       pars[assign == i], x, ...))))
      })

