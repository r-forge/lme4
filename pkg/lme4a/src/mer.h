// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "MatrixNs.h"

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

	MatrixNs::chmFa L;
	MatrixNs::chmSp Lambda, Ut, Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta, u;
	double *ldL2;
    };
    
    class merResp {
    public:
	merResp(Rcpp::S4);
	void updateL(reModule&);
	//<  cu <- solve(L, solve(L, crossprod(Lambda, Utr), sys = "P"), sys = "L")
	void updateWrss();	//< resid := y - mu; wrss := sum((sqrtrwts * resid)^2)
	
	Rcpp::NumericVector Utr, Vtr, cbeta,
	    cu, mu, offset, resid, weights, y; 
	MatrixNs::chmDn cUtr, ccu;
	double *wrss;
    };

    class feModule {
    public:
	feModule(Rcpp::S4 xp) :
	    beta(SEXP(xp.slot("beta"))) {
	    Rcpp::NumericVector ldR2Vec(SEXP(xp.slot("ldR2")));
	    ldR2 = ldR2Vec.begin();
	}

	double *ldR2;
	Rcpp::NumericVector beta;
    };
    
    class deFeMod : public feModule {
    public:
	deFeMod(Rcpp::S4 xp) :
	    feModule(xp),
	    X(Rcpp::S4(SEXP(xp.slot("X")))),
	    RZX(Rcpp::S4(SEXP(xp.slot("RZX")))),
	    RX(Rcpp::S4(SEXP(xp.slot("RX")))),
	    cRZX(RZX) { }
	
	MatrixNs::dgeMatrix X, RZX;
	MatrixNs::Cholesky RX;
	MatrixNs::chmDn cRZX;
    };

    class spFeMod : public feModule {
    public:
	spFeMod(Rcpp::S4 xp) :
	    feModule(xp),
	    X(Rcpp::S4(SEXP(xp.slot("X")))),
	    RZX(Rcpp::S4(SEXP(xp.slot("RZX")))),
	    RX(Rcpp::S4(SEXP(xp.slot("RX")))) { }
	
	MatrixNs::chmSp X, RZX;
	MatrixNs::chmFa RX;
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
	
	MatrixNs::chmSp ZtX;
	MatrixNs::chmFa XtX;
    };

    class lmerDeFeMod : public deFeMod {
    public:
	lmerDeFeMod(Rcpp::S4 xp) :
	    deFeMod(xp),
	    ZtX(Rcpp::S4(SEXP(xp.slot("ZtX")))), cZtX(ZtX),
	    XtX(Rcpp::S4(SEXP(xp.slot("XtX")))) { }

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
	MatrixNs::chmDn cZtX;
	double *ldR2;
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
	//< ldL2 + ldR2 + (n - p) * (1 + log(2 * pi * pwrss/(n - p)))
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
	//< ldL2 + ldR2 + (n - p) * (1 + log(2 * pi * pwrss/(n - p)))
	double updateTheta(const Rcpp::NumericVector&);
	lmerDeFeMod fe;
    };
}

#endif

