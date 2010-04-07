stopifnot(require(lme4a),
          data(ergoStool, package="MEMSS") == "ergoStool")

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

