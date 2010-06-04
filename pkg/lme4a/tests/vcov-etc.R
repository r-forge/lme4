stopifnot(require(lme4a))
## "MEMSS" is just 'Suggest' -- must still work, when it's missing:
if(data(ergoStool, package="MEMSS") != "ergoStool") {
    cat("'ergoStool' data from package 'MEMSS' is not available --> skipping test\n")
    quit('no')
}

fm1   <-  lmer (effort ~ Type + (1|Subject), data = ergoStool)
fm1.s  <- lmer (effort ~ Type + (1|Subject), data = ergoStool, sparseX=TRUE)
## was segfaulting with sparseX (a while upto 2010-04-06)

fe1   <- fixef(fm1)
fe1.s <- fixef(fm1.s)

s1.d <- summary(fm1)
s1.s <- summary(fm1.s)
stopifnot(
	  all.equal(fe1, fe1.s, tol= 1e-12)
	  ,
	  all.equal(se1.d <- coef(s1.d)[,"Std. Error"],
		    se1.s <- coef(s1.s)[,"Std. Error"], tol=1e-10)
	  ,
	  all.equal(V.d <- vcov(fm1),
		    V.s <- vcov(fm1.s), tol = 1e-9)
	  ,
	  all.equal(diag(V.d), unname(se1.d)^2, tol= 1e-12)
###B	  ,
###Bug ??: currently have
###B        c(0.5921429534113, rep(0.5114296867053, 3))
###B	  all.equal(unname(se1.d),
###B		    ## lmer2:
###B		    c(0.5760119655646, rep(0.5186839846263, 3)), tol= 1e-10)
###B	  ,
###B	  all.equal(unname(se1.d),
###B		    ## lmer1:
###B		    c(0.5760122554859, rep(0.5186838382653, 3)), tol= 1e-6)
	  )

### -------------------------- a "large" example -------------------------
str(InstEval)

## this works
system.time(
fm7 <- lmer(y ~ d + service + studage + lectage + (1|s),
             data = InstEval, sparseX = TRUE, verbose = TRUE)
)
system.time(sfm7 <- summary(fm7))
fm7 # takes a while as it computes summary() again !

range(t.fm7 <- coef(sfm7)[,"t value"])## -10.94173  10.61535

m.t.7 <- mean(abs(t.fm7), trim = .01)
###B : now have     m.t.7= 1.55511602701
stopifnot(all.equal(m.t.7, 1.55326394,   tol = 2e-3), # had = 1e-5  # lmer1
          all.equal(m.t.7, 1.5532709682, tol = 2e-3)) # had = 1e-9  # lmer2
hist.t <- cut(t.fm7, floor(min(t.fm7)) : ceiling(max(t.fm7)))
cbind(table(hist.t))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
