library(Gqr)
wfmt <- "http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_hermite_probabilist/hermite_probabilist_%03d_w.txt "
xfmt <- "http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_hermite_probabilist/hermite_probabilist_%03d_x.txt "
nn <- 2^(1:7) - 1L
for (n in nn) {
    loc <- GaussQuad(n, "H", b=0.5)
    webw <- scan(url(sprintf(wfmt, n)))
    webx <- scan(url(sprintf(xfmt, n)))
    stopifnot(all.equal(webw, loc$weights), all.equal(webx, loc$knots))
}
