GaussQuad <-
    function(n, rule = c("Legendre", "Chebyshev", "Gegenbauer", "Jacobi",
                         "Laguerre", "Hermite", "Exponential", "Rational"),
             a=0, b=1, alpha=0, beta=0)
{
    ## Note: rule "Type2" is not yet supported in C code
    rules <- eval(formals()$rule)
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

integrateGQ <- function(f, lower, upper, ...,
                        n = 13, rule = "Legendre",
                        a=lower, b=upper, alpha=0, beta=0)
{
    f <- match.fun(f)
    stopifnot(length(a) == 1, is.finite(a),
              length(b) == 1, is.finite(b),
              length(alpha) == 1, is.finite(alpha),
              length(beta) == 1, is.finite(beta))
    rules <- eval(formals(GaussQuad)$rule)
    rule <- match.arg(rule, choices = rules)
    GQ <- GaussQuad(n, rule=rule, a=a, b=b, alpha=alpha, beta=beta)
    with(GQ, sum(weights * f(knots, ...)))
}
