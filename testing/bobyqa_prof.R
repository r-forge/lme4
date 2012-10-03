prof2 <- function(theta) {
  pp <- get("pp",environment(ff))
  ff2 <- function(beta) {
    ff(c(theta,beta))
  }
  optval <- optimx(pp$beta(0),
                   ff2, method="bobyqa", control=list(rhobeg=2e-4,rhoend=2e-7))
  optval$fvalues[[1]]
}

prof2(val[1])

fn <- "bobyqa1_prof.RData"
vvec <- seq(0,1.2,by=0.025)
pp <- numeric(length(vvec))
for (i in seq_along(pp)) {
  cat(i,vvec[i],"\n")
  pp[i] <- prof2(vvec[i])
}
save("pp","vvec",file=fn)


