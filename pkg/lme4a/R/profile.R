prTheta <- function(fenv, thmin = 0, thmax = 5, n = 101)
{
    thmin <- as.numeric(thmin)[1]
    thmax <- as.numeric(thmax)[1]
    n <- as.integer(n)[1]
    paropt <- getPars(fenv)
    stopifnot(is.environment(fenv),
              length(paropt) == 1,
              0 <= thmin, thmin < thmax,
              n > 1)
    thseq <- seq(thmin, thmax, len.out = n)
    data.frame(theta = thseq, dev = sapply(thseq, setPars, x = fenv))
}

