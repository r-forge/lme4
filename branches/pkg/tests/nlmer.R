library(lme4)

showProc.time <- function() { ## CPU elapsed __since last called__
    .ot <- .pc
    .pc <<- proc.time()
    cat('Time elapsed: ', (.pc - .ot)[1:3],'\n')
}
allEQ <- function(x,y, tolerance = 4e-4, ...)
    all.equal.numeric(x,y, tolerance=tolerance, ...)
.pc <- proc.time()

## 'Theoph' Data modeling

Th.start <- c(lKe = -2.5, lKa = 0.5, lCl = -3)
(nm2 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKe+lKa+lCl|Subject),
              Theoph, start = Th.start))
showProc.time() # ~ 5.7s {dual-opteron 2814, on 64b, no optim.}

(nm3 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKe|Subject) + (lKa|Subject) + (lCl|Subject),
              Theoph, start = Th.start))
showProc.time() # ~ 3.2s

## dropping   lKe  from random effects:
(nm4 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKa+lCl|Subject),
              Theoph, start = Th.start))
showProc.time()

(nm5 <- nlmer(conc ~ SSfol(Dose, Time,lKe, lKa, lCl) ~
              (lKa|Subject) + (lCl|Subject),
              Theoph, start = Th.start))
showProc.time()

e3 <- expand(nm3)
stopifnot(identical(sapply(e3, class),
                    c(sigma = "numeric", P = "pMatrix",
                      T = "dtCMatrix", S = "ddiMatrix"))
          , allEQ(e3$sigma, c(sigmaML = 0.76921295))
          , all(e3$P@perm == outer(12*(0:2), 1:12, "+"))
          , identical(as(e3$T, "diagonalMatrix"), Diagonal(3*12))
          , allEQ(e3$S@x, rep(c(0, 0.620071, 0.163092), each=12))
          )

e2 <- expand(nm2) # -> gave error!
stopifnot(identical(sapply(e2, class),
                    c(sigma = "numeric", P = "pMatrix",
                      T = "dtCMatrix", S = "ddiMatrix"))
          , allEQ(e2$sigma, c(sigmaML = 0.7692315))
          , all(e2$P@perm == outer(12*(0:2), 1:12, "+"))
          , diag(e2$T == 1), nnzero(e2$T) == 36 + 24 + 12
          , allEQ(unique(e2$T@x),
                  c(1, 0.001, 0.0564, -0.00328), tol = .009)
          , allEQ(e2$S@x, rep(c(1e-6, 0.6198669, 0.1629919), each=12))
          )

showProc.time()


