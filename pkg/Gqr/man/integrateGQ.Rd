\name{integrateGQ}
\alias{integrateGQ}
\title{Numerical Integration via Gauss-Hermite}
\description{
  Numerical integration of one-dimensional functions, using
  Gauss-Hermite, with respect to eight different kernel (\code{rule}s).
}
\usage{
integrateGQ(f, lower, upper, ..., n = 13, rule = "Legendre",
            a = lower, b = upper, alpha = 0, beta = 0)
}
\arguments{
  \item{f}{an \R \code{\link{function}} taking a numeric first argument,
    the \sQuote{integrand} of the numerical integration (aka
    \dQuote{quadrature}).}
  \item{lower, upper}{the limits of integration.  Currently \emph{must}
    be finite.}
  \item{\dots}{additional arguments to be passed to \code{f}.}

  \item{n}{number of Gauss-Hermite quadrature points, see
    \code{\link{GaussQuad}}.}
  \item{rule}{a character string specifying one of the quadrature rules,
    see \code{\link{GaussQuad}}.}
  \item{a,b, alpha, beta}{numeric scalars, specifying the integral to be
    computed; the meaning of these depends much on the \code{rule}, see
    the \code{Details} in  \code{\link{GaussQuad}}.}
}
\value{
  a number, simply
  \deqn{sum_{i=1}^n w_i f(x_i),}{sum(i=1:n; w[i] * f(x[i])),}
  where \eqn{x_i}{x[i]} are the \code{knots}
  and   \eqn{w_i}{w[i]} are the \code{weights} returned from
  \code{\link{GaussQuad}(n, rule, ..)}.
}
%%\author{Martin Maechler}
\seealso{
  \code{\link{integrate}()}, \R's \dQuote{built-in} numerical integration.
  \code{\link{GaussQuad}}, also for references.
}
\examples{
 f1 <- function(x) x^2
 (i1 <- integrate  (f1, -1, 3))
 (i2 <- integrateGQ(f1, -1, 3))
 stopifnot(all.equal(i2, 28/3, tol = 1e-12), all.equal(i1$value, i2))
 ## n does not matter for such a simple f(): always ~ 14 digits accuracy:
 for(n in 2:20) cat(sprintf("\%3d: \%5.2f\n", n,
                    log10(abs(28/3 - integrateGQ(f1, -1, 3, n=n)))))

 ## However, it *does* for Hermite, and a "hard" case:
 gi3 <- function(x) 1/(1+ abs(x)^3)
 n. <- c(3:130,150, 170, 200, 250, 350, 500, 1000)
 Eg. <- sapply(n., function(n)
               integrateGQ(gi3, n=n, rule = "Hermite", a=0, b=1/2))
 ## (it *is* very fast)
 plot(Eg. ~ n., type = "l", col=2,
      subset = n. <= 130, ylim = c(1.65, 1.73))

 ## "Classical" numerical integration:
 (Icl <- integrate(function(u) gi3(u)*exp(-u^2/2), -Inf,Inf,
                   rel.tol= 1e-13))
 ## 1.6810227521877 with absolute error < 5.2e-14
 Ic <- Icl$value #  (let's assume it's "true")

 xax <- require("sfsmisc")
 plot(abs(Eg. - Ic) ~ n., type = "o", col=2, log="xy",
      axes = !xax,
      main = "Gauss-Hermite approximation error [log-log scale]")
 if(xax) { eaxis(1); eaxis(2) }

 ## Compute  E[ g(X) ]  where  X ~ N(0, 1)  via  Gauss-Hermite
 E.norm <- function(f, n = 32)
    integrateGQ(f, n=n, rule = "Hermite", a=0, b=1/2) / sqrt(2*pi)

 integrands <- list(
    g1 = function(x) (x-2)^2,
    g2 = function(x) (x+3)^3,
    g4 = function(x) x^4,
    gg = function(x) 1/(1+ abs(x)^3) )

 op <- par(mfrow= c(2,2), mgp=c(1.5,.6, 0), mar=.1+c(3,3,3,0))
 for(i in seq_along(integrands)) {
    g <- integrands[[i]]
    nm <- names(integrands)[i]
    cat(nm,":", deparse(expr <- body(g)), ":\n")
    curve(g, -3, 3, main = expr); abline(h=0, lty=3)
    cat("classical integrate():\n")
    print(i1 <- integrate(function(u) g(u)*dnorm(u), -Inf, Inf,
                          rel.tol = 1e-13),  digits=12)
    cat("\n Gauss-Hermite (n= 5):", i2 <- E.norm(g, n= 5),"\n")
    cat(  " Gauss-Hermite (n=16):", i3 <- E.norm(g, n=16),"\n")
    cat(  " Gauss-Hermite (n=64):", i4 <- E.norm(g, n=64),"\n")
    cat(  " Gauss-Hermite (n=128):",i5 <- E.norm(g,n=128),"\n\n")
    legend("top",
           c(paste("E[g(X)] {integr.()} = ", formatC(i1$value)),
             paste("E[g(X)] {Hermite_5 } =", formatC(i2)),
             paste("E[g(X)] {Hermite_16} =", formatC(i3)),
             paste("E[g(X)] {Hermite_64} =", formatC(i4)),
             paste("E[g(X)] {Hermite_128}=", formatC(i5))),
           bty="n")
 }
 par(op)

 doExtras <- interactive()## takes a bit time: using *large* n:
 if(doExtras) {
   gg <- integrands$gg
   nn <- 2^seq(2,13, by=.5)
   E.gg <- sapply(nn, E.norm, f = gg)
   do.eax <- require("sfsmisc")
   ylb <- expression(E[N(0,1)](g(X)))
   plot(E.gg ~ nn, log = "x", col=2, type="b", xaxt= if(do.eax) "n",
        main = substitute(list(EGX, ~~g(x) == GG),
                          list(EGX = ylb[[1]], GG = body(gg))),
        ylab = "E[g(X)]")
   if(do.eax) eaxis(1)
   print(cbind(round(nn,2), E.gg), digits = 12)

   plot(E.gg - i1$val ~ nn, subset = nn >= 100, log="x", type="b", col=2)
   abline(h=0,lty=3)
   plot(abs(E.gg - i1$val) ~ nn, subset = nn >= 100, log="xy", type="b", col=2)
 }
}
\keyword{math}
\keyword{utilities}
