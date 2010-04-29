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
	  ,
	  all.equal(unname(se1.d),
		    c(0.5760122554859, rep(0.5186838382653, 3)), tol= 1e-10)
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

stopifnot(all.equal(mean(abs(t.fm7), trim = .01), 1.55326394))

hist.t <- cut(t.fm7, floor(min(t.fm7)) : ceiling(max(t.fm7)))
cbind(table(hist.t))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
