// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "MatrixNs.h"
#include "glmFamily.h"

namespace mer {
    class merResp;		// forward declaration

    class reModule {
    public:
	reModule(Rcpp::S4);

	void updateTheta(const Rcpp::NumericVector&);
	//< Lambda@x[] <- theta[Lind]; Ut <- crossprod(Lambda, Zt); update(L,Ut,1); ldL2
	double sqLenU() const {	//< squared length of u
	    return std::inner_product(u.begin(), u.end(), u.begin(), double());
	}
	void updateU(const merResp&);
	//< u <- solve(L, solve(L, cu - RZX %*% beta, sys = "Lt"), sys = "Pt")
	void incGamma(Rcpp::NumericVector&);
	//< gamma += crossprod(Ut, u)

	MatrixNs::chmFr L;
	MatrixNs::chmSp Lambda, Ut, Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta, u;
	double *ldL2;
    };

    class rwReMod : public reModule {
    public:
	rwReMod(Rcpp::S4 xp) :
	    reModule(xp),
	    ubase(SEXP(xp.slot("ubase"))) {}
	void incGamma(Rcpp::NumericVector&);
	//< gamma += crossprod(Ut, u + ubase)
	void updateUt(const merResp&);

	Rcpp::NumericVector ubase;
    };

    class merResp {
    public:
	merResp(Rcpp::S4);
	void updateL(reModule&);
	//<  cu <- solve(L, solve(L, crossprod(Lambda, Utr), sys = "P"), sys = "L")
	double updateWrss(); //< wtres <- sqrtrwts * (y - mu); wrss <- sum(wtres^2)

	Rcpp::NumericVector Utr, Vtr, cbeta,
	    cu, mu, offset, wtres, weights, y;
	MatrixNs::chmDn cUtr, ccu;
	double *wrss;
    };

    class feModule {
    public:
	feModule(Rcpp::S4 xp) :
	    beta(SEXP(xp.slot("beta"))) {
	    Rcpp::NumericVector l_vec(SEXP(xp.slot("ldRX2")));
	    ldRX2 = l_vec.begin();
	}

	double *ldRX2;
	Rcpp::NumericVector beta;
    };

    class deFeMod : public feModule {
    public:
	deFeMod(Rcpp::S4 xp) :
	    feModule(xp)
	    , X(Rcpp::S4(SEXP(xp.slot("X"))))
	    , RZX(Rcpp::S4(SEXP(xp.slot("RZX"))))
	    , RX(Rcpp::S4(SEXP(xp.slot("RX"))))
	{ }

	MatrixNs::dgeMatrix X, RZX;
	MatrixNs::Cholesky RX;
    };

    class lmerDeFeMod : public deFeMod {
    public:
	lmerDeFeMod(Rcpp::S4 xp) :
	    deFeMod(xp)
	    , ZtX(Rcpp::S4(SEXP(xp.slot("ZtX"))))
	    , XtX(Rcpp::S4(SEXP(xp.slot("XtX"))))
	{ }

	void updateRzxRx(reModule&);
	//< RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), sys = "P"), sys = "L")
        //< RX <- chol(XtX - crossprod(RZX))
	void updateBeta(merResp&);
	//< resp@cbeta <- Vtr - crossprod(RZX, cu)
	//< beta <- solve(RX, solve(t(RX), resp@cbeta))
	//< resp@cu <- resp@cu - RZX %*% beta
	void incGamma(Rcpp::NumericVector &gam) {
	    X.dgemv('N', 1., beta, 1., gam);
	}
	//< gamma += crossprod(Ut, u)

	MatrixNs::dgeMatrix ZtX;
	MatrixNs::dpoMatrix XtX;
    };

    class spFeMod : public feModule {
    public:
	spFeMod(Rcpp::S4 xp) :
	    feModule(xp),
	    X(Rcpp::S4(SEXP(xp.slot("X")))),
	    RZX(Rcpp::S4(SEXP(xp.slot("RZX")))),
	    RX(Rcpp::S4(SEXP(xp.slot("RX")))) { }

	MatrixNs::chmSp X, RZX;
	MatrixNs::chmFr RX;
    };

    class lmerSpFeMod : public spFeMod {
    public:
	lmerSpFeMod(Rcpp::S4 xp) :
	    spFeMod(xp),
	    ZtX(Rcpp::S4(SEXP(xp.slot("ZtX")))),
	    XtX(Rcpp::S4(SEXP(xp.slot("XtX")))) { }

	void updateRzxRx(reModule&);
	//< RZX <<- solve(L, solve(L, crossprod(Lambda, ZtX), sys = "P"), sys = "L")
        //< RX <- chol(XtX - crossprod(RZX))
	void updateBeta(merResp&);
	//< resp@cbeta <- Vtr - crossprod(RZX, cu)
	//< beta <- solve(RX, solve(t(RX), resp@cbeta))
	//< resp@cu <- resp@cu - RZX %*% beta
	void incGamma(Rcpp::NumericVector &gam) {
	    MatrixNs::chmDn bb(beta), gg(gam);
	    X.dmult('N', 1., 1., bb, gg);
	}
	//< gamma += crossprod(Ut, u)

	MatrixNs::chmSp ZtX, XtX;
    };

    class lmer {
    public:
	lmer(Rcpp::S4 xp) :
	    re(Rcpp::S4(SEXP(xp.slot("re")))),
	    resp(Rcpp::S4(SEXP(xp.slot("resp")))),
	    REML(SEXP(xp.slot("REML"))) {
	    reml = (bool)*REML.begin();
	}

	double deviance();
	//< ldL2 + n *(1 + log(2 * pi * pwrss/n))

	reModule re;
	merResp resp;
	Rcpp::LogicalVector REML;
	bool reml;
    };

    class lmerDe : public lmer {
    public:
	lmerDe(Rcpp::S4 xp) :
	    lmer(xp),
	    fe(Rcpp::S4(SEXP(xp.slot("fe")))) {
	}
	double reCrit();
	//< ldL2 + ldRX2 + (n - p) * (1 + log(2 * pi * pwrss/(n - p)))
	double updateTheta(const Rcpp::NumericVector&);
	lmerDeFeMod fe;
    };

    class lmerSp : public lmer {
    public:
	lmerSp(Rcpp::S4 xp) :
	    lmer(xp),
	    fe(Rcpp::S4(SEXP(xp.slot("fe")))) {
	}
	double reCrit();
	//< ldL2 + ldRX2 + (n - p) * (1 + log(2 * pi * pwrss/(n - p)))
	double updateTheta(const Rcpp::NumericVector&);
	lmerSpFeMod fe;
    };

    class rwResp : public merResp {
    public:
	rwResp(Rcpp::S4 xp) :
	    merResp(xp),
	    gamma(SEXP(xp.slot("gamma"))),
	    sqrtrwt(SEXP(xp.slot("sqrtrwt"))),
	    sqrtXwt(Rcpp::S4(SEXP(xp.slot("sqrtXwt")))) {}
	double updateWrss();
	Rcpp::NumericVector gamma, sqrtrwt;
	MatrixNs::dgeMatrix sqrtXwt;
    };

    class glmerResp : public rwResp {
    public:
	glmerResp(Rcpp::S4 xp) :
	    rwResp(xp),
	    family(SEXP(xp.slot("family"))),
	    muEta(SEXP(xp.slot("muEta"))),
	    n(SEXP(xp.slot("n"))),
	    var(SEXP(xp.slot("var"))) {}

	void linkFun(){family.linkFun(gamma, mu);}
	void linkInv(){family.linkInv(mu, gamma);}
	void MuEta(){family.muEta(muEta, gamma);}
	void variance(){family.variance(var, mu);}
	void updateSqrtRWt();
	void updateSqrtXWt();

	glmFamily family;
	Rcpp::NumericVector muEta, n, var;
    };
    
    class rwDeFeMod : public deFeMod {
    public:
	rwDeFeMod(Rcpp::S4 xp) :
	    deFeMod(xp),
	    V(Rcpp::S4(SEXP(xp.slot("V")))),
	    betabase(SEXP(xp.slot("betabase"))) {}

	void incGamma(Rcpp::NumericVector &gam);
//	void updateV(rwResp const&);
	void updateV(rwResp&);
	MatrixNs::dgeMatrix V;
	Rcpp::NumericVector betabase;
    };

    class rwSpFeMod : public spFeMod {
    public:
	rwSpFeMod(Rcpp::S4 xp) :
	    spFeMod(xp),
	    V(Rcpp::S4(SEXP(xp.slot("V")))),
	    betabase(SEXP(xp.slot("betabase"))) {}

	MatrixNs::chmSp V;
	Rcpp::NumericVector betabase;
    };
    
    class glmer {
    public:
	glmer(Rcpp::S4 xp) :
	    re(Rcpp::S4(SEXP(xp.slot("re")))),
	    resp(Rcpp::S4(SEXP(xp.slot("resp")))) {}
	
	rwReMod re;
	glmerResp resp;
    };
    
    class glmerDe : public glmer {
    public:
	glmerDe(Rcpp::S4 xp) :
	    glmer(xp),
	    fe(Rcpp::S4(SEXP(xp.slot("fe")))) {}

	void updateGamma();
	double IRLS();

	rwDeFeMod fe;
    };

    class glmerSp : public glmer {
    public:
	glmerSp(Rcpp::S4 xp) :
	    glmer(xp),
	    fe(Rcpp::S4(SEXP(xp.slot("fe")))) {}

	rwSpFeMod fe;
    };

}

#endif







