sparseRasch <- function(dat, maxirls = 200, tol = 1e-5)
{
    ## Convert individual columns of a data frame
    as.binaryNumeric <- function(x) {
        if (is.logical(x)) return(as.numeric(x))
        if (is.factor(x)) {
            if (length(levels(x)) != 2)
                stop("factors in the data frame must have exactly 2 levels")
            return(as.numeric(x) - 1L)
        }
    }
    
    ## Massage the data into a numeric, binary matrix
    if (is.data.frame(dat))
        dat <- do.call("cbind", lapply(dat, as.binaryNumeric))
    dat <- as.matrix(dat)
    storage.mode(dat) <- "double"
    nr <- nrow(dat)
    nc <- ncol(dat)
    y <- as.vector(dat)
    stopifnot(length(unique(y)) == 2, min(y) == 0, max(y) == 1)

    ## Generate the transpose of the model matrix using the
    ## contr.treatment scheme.  Ability parameters come first then
    ## difficulties then the intercept, which is the logit of the
    ## probability of a correct response by the first subject on the
    ## first question.
    MM <- rbind2(rbind2(as(gl(nr, 1, nr * nc), "sparseMatrix")[-1, ],
                        as(gl(nc, nr), "sparseMatrix")[-1, ]),
                 t(rep(1, nr * nc)))
    M1 <- MM
    
    beta <- rnorm(nrow(MM), mean = 0, sd = 0.1)
    for (i in c(0, seq_len(maxirls))) { # IRLS iterations
        eta <- crossprod(MM, beta)@x
        beta_old <- beta
        mu <- .Call("logit_linkinv", eta, PACKAGE = "stats")
        swts <- sqrt(1/(mu * (1 - mu)))
        wtres <- swts * (y - mu)
        M1@x <- swts * .Call("logit_mu_eta", eta, PACKAGE = "stats")
        L <- Cholesky(tcrossprod(M1), perm = FALSE, LDL = FALSE)
        inc1 <- solve(L, wtres, "L")
        crit <- sqrt((inc1 %*% inc1)/(wtres %*% wtres))
        beta <- beta + solve(L, inc1, "Lt")
        if (crit < tol) break
    }
    MM@x[] <- 1
    list(beta = beta, Xt = MM, mu = mu, L = L)
}
