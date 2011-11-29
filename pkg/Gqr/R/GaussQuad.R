GaussQuad <-
    function(n, rule, a=0, b=1, alpha=0, beta=0) {
        rules <- c("Legendre", "Chebyshev", "Gegenbauer", "Jacobi",
                   "Laguerre", "Hermite", "Exponential", "Rational", "Type2")
        if (is.na(rr <- pmatch(rule, rules)))
            stop("Unknown quadrature rule, ", rule)
        stopifnot(length(n <- as.integer(n)) == 1L,
                  length(a <- as.numeric(a)) == 1L,
                  length(b <- as.numeric(b)) == 1L,
                  length(alpha <- as.numeric(alpha)) == 1L,
                  length(beta <- as.numeric(beta)) == 1L,
                  n > 0L)
        .Call(Gauss_Quad, n, rr, a, b, alpha, beta)
    }

