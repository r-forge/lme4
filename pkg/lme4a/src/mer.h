// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "MatrixNs.h"
#include "glmFamily.h"

namespace mer {
    class merResp;		// forward declaration

    class reModule {
	MatrixNs::chmFr L;
	MatrixNs::chmSp Lambda, Ut, Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta;
	double *d_ldL2;
	Rcpp::NumericVector d_u;
    public:
	reModule(Rcpp::S4);

	CHM_SP SupdateL(MatrixNs::chmSp const&) const;
	void setU(const double *);
	const Rcpp::NumericVector &u() const;
	double sqLenU() const;
	double ldL2() const;
	void DupdateL(MatrixNs::chmDn const&,MatrixNs::chmDn&) const;
	void incGamma(Rcpp::NumericVector&) const;
	void rwUpdateL(MatrixNs::dgeMatrix const&,
		       Rcpp::NumericVector const&,
		       Rcpp::NumericVector const&,
		       double*);
	void updateTheta(Rcpp::NumericVector const&);
	void updateU(merResp const&);
    };

    class merResp {
    public:
	merResp(Rcpp::S4);
	void updateL(reModule const&);
	void updateMu(Rcpp::NumericVector const&);
	double updateWrss();

	Rcpp::NumericVector Utr, Vtr, cbeta, cu, mu,
	    offset, wtres, weights, y;
	double *wrss;
    };

    class feModule {
    protected:
	Rcpp::NumericVector d_beta, d_ldRX2;
    public:
	feModule(Rcpp::S4 xp);
	void setBeta(const double*);
	const Rcpp::NumericVector& beta() const;
	double ldRX2() const;
    };

    class deFeMod : public feModule {
    protected:
	MatrixNs::dgeMatrix X, RZX;
	MatrixNs::Cholesky d_RX;
    public:
	deFeMod(Rcpp::S4 xp);
	const MatrixNs::Cholesky &RX() const;
	void incGamma(Rcpp::NumericVector &gam) const;
    };

    class lmerDeFeMod : public deFeMod {
	MatrixNs::dgeMatrix ZtX;
	MatrixNs::dpoMatrix XtX;
    public:
	lmerDeFeMod(Rcpp::S4 xp);

	void updateRzxRx(reModule const&);
	void updateBeta(merResp&);
    };

    class spFeMod : public feModule {
    protected:
	MatrixNs::chmSp X, RZX;
	MatrixNs::chmFr RX;
    public:
	spFeMod(Rcpp::S4 xp);
	void incGamma(Rcpp::NumericVector &gam) const;
    };

    class lmerSpFeMod : public spFeMod {
	MatrixNs::chmSp ZtX, XtX;
    public:
	lmerSpFeMod(Rcpp::S4 xp);
	void updateRzxRx(reModule const&);
	void updateBeta(merResp&);
    };

    class rwResp : public merResp {
    public:
	rwResp(Rcpp::S4 xp);
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
	rwDeFeMod(Rcpp::S4 xp);
	void updateV(rwResp const&);
	void updateRX(bool);
	void dpotrs(Rcpp::NumericVector&) const;
	MatrixNs::dgeMatrix V;
    };

    class rwSpFeMod : public spFeMod {
    public:
	rwSpFeMod(Rcpp::S4 xp) :
	    spFeMod(xp),
	    V(Rcpp::S4(xp.slot("V"))){};
	MatrixNs::chmSp V;
    };
    
    class glmer {
    protected:
	reModule re;
	glmerResp resp;
	int verbose;
    public:
	glmer(Rcpp::S4 xp);
	double Laplace();
    };
    
    class glmerDe : public glmer {
	rwDeFeMod fe;
    public:
	glmerDe(Rcpp::S4 xp);
	void updateGamma();
	double IRLS(int);
	double PIRLS(int);
    };

    class glmerSp : public glmer {
	rwSpFeMod fe;
    public:
	glmerSp(Rcpp::S4 xp);
	void updateGamma();
	double IRLS(int);
    };

    template<typename T>
    class lmer {
	reModule re;
	T fe;
	merResp resp;
	bool reml;
	static double l2PI;
    public:
	lmer(Rcpp::S4 xp);

	double deviance();
	double reCrit();
	double updateTheta(const Rcpp::NumericVector&);
    };
    
    template<typename T>
    inline lmer<T>::lmer(Rcpp::S4 xp) :
	re(Rcpp::S4(SEXP(xp.slot("re")))),
	fe(Rcpp::S4(xp.slot("fe"))),
	resp(Rcpp::S4(SEXP(xp.slot("resp"))))
    {
	Rcpp::LogicalVector REML = xp.slot("REML");
	reml = (bool)*REML.begin();
    }

    template<typename T>
    double lmer<T>::l2PI = log(2. * PI);
    
/** 
 * Evaluate the deviance
 * 
 * @return ldL2 + n *(1 + log(2 * pi * pwrss/n))
 */
    template<typename T> inline
    double lmer<T>::deviance() {
	double nn = (double)resp.y.size(),
	    prss = re.sqLenU() + *resp.wrss;
	return re.ldL2() + nn * (1 + l2PI + log(prss/nn));
    }

/** 
 * Evaluate the REML criterion
 * 
 * @return ldL2 + ldRX2 + (n - p) * (1 + log(2 * pi * pwrss/(n - p)))
 */
    template<typename T> inline
    double lmer<T>::reCrit() {
	double nmp = (double)(resp.y.size() - fe.beta().size()),
	    prss = re.sqLenU() + *resp.wrss;
	return re.ldL2()+fe.ldRX2()+nmp*(1 + l2PI + log(prss/nmp));
    }

    template<typename T> inline
    double lmer<T>::updateTheta(Rcpp::NumericVector const &nt) {
	re.updateTheta(nt);
	resp.updateL(re);
	fe.updateRzxRx(re);
	fe.updateBeta(resp);
	re.updateU(resp);

	Rcpp::NumericVector gamma(resp.offset.size());
	std::copy(resp.offset.begin(), resp.offset.end(), gamma.begin());
	fe.incGamma(gamma);
	re.incGamma(gamma);
	resp.updateMu(gamma);
	
	resp.updateWrss();	// update resp.wtres and resp.wrss
	return reml ? reCrit() : deviance();
    }
}

#endif
