library(Gqr)

## Do not do the checks per default "on CRAN"
(doExtras <- {
    interactive() ||
        grepl("bates", Sys.getenv("USER")) ||
        nzchar(Sys.getenv("R_Gqr_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
})

if(doExtras) {
    urldir <- "http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_hermite_probabilist/"
    wfmt <- paste0(urldir,"hermite_probabilist_%03d_w.txt")
    xfmt <- paste0(urldir,"hermite_probabilist_%03d_x.txt")
}
nn <- 2^(1:7) - 1L
for (n in nn) {
    loc <- GaussQuad(n, "H", b=0.5)
    if(doExtras) {
        webw <- scan(url(sprintf(wfmt, n)))
        webx <- scan(url(sprintf(xfmt, n)))
        stopifnot(all.equal(webw, loc$weights), all.equal(webx, loc$knots))
    }
}
