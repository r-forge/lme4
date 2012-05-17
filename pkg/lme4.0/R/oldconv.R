convert_old_lme4 <- function(x) {
    cc <- class(x)
    attr(cc,"package") <- "lme4.0"
    class(x) <- cc
    x
}
