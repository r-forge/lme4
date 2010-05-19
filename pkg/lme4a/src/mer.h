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
	void incGamma(double*) const;
	void rwUpdateL(Rcpp::NumericMatrix const&,
		       Rcpp::NumericVector const&,
		       Rcpp::NumericVector const&,
		       double*);
	void updateTheta(Rcpp::NumericVector const&);
	void updateU(merResp const&);
    };

    class merResp {
    protected:
	Rcpp::NumericVector d_sqrtrwt;
    public:
	merResp(Rcpp::S4);
	void updateL(reModule const&);
	void updateMu(Rcpp::NumericVector const&);
	double updateWrss();
	const double *sqrtrwt() const;

	Rcpp::NumericVector Utr, Vtr, cbeta, cu, mu,
	    offset, wtres, weights, y;
	double *wrss;
    };

    class feModule {
    protected:
	Rcpp::NumericVector d_beta, d_ldRX2;
    public:
	feModule(Rcpp::S4 xp);
	void setBeta(Rcpp::NumericVector const&);
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
	void incGamma(double*) const;
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
    protected:
	Rcpp::NumericVector d_gamma;
	Rcpp::NumericMatrix d_sqrtXwt;
    public:
	rwResp(Rcpp::S4 xp);
	const Rcpp::NumericMatrix &sqrtXwt() const;
	double *gamma();
    };

    class glmerResp : public rwResp {
    protected:
	glmFamily family;
	Rcpp::NumericVector muEta, n, d_var;
    public:
	glmerResp(Rcpp::S4 xp);
	const Rcpp::NumericVector &var() const {return d_var;}
	void linkFun(){family.linkFun(d_gamma, mu);}
	void linkInv(){family.linkInv(mu, d_gamma);}
	void MuEta(){family.muEta(muEta, d_gamma);}
	void variance(){family.variance(d_var, mu);}
	Rcpp::NumericVector devResid(){return family.devResid(mu, weights, y);}
	void updateSqrtRWt();
	void updateSqrtXWt();
    };
    
    class nlmerResp : public rwResp {
	Rcpp::NumericVector d_eta;
	Rcpp::Environment nlenv;
	Rcpp::Language nlmod;
	Rcpp::CharacterVector pnames;
    public:
	nlmerResp(Rcpp::S4 xp);
	const Rcpp::NumericVector &eta() const;
	void nlEval();
//	const Rcpp::NumericVector &sqrtXwt() const;
    };
    
    class rwDeFeMod : public deFeMod {
	MatrixNs::dgeMatrix d_V;
    public:
	rwDeFeMod(Rcpp::S4 xp);
	void updateV(Rcpp::NumericMatrix const&);
	void updateRX(bool);
	void dpotrs(Rcpp::NumericVector&) const;
	const MatrixNs::dgeMatrix &V() const;
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

    class nlmer {
    protected:
	reModule re;
	nlmerResp resp;
    public:
	nlmer(Rcpp::S4 xp);
	void nlEval();
	double Laplace();
    };
    
    class nlmerDe : public nlmer {
	rwDeFeMod fe;
    public:
	nlmerDe(Rcpp::S4 xp);
	void updateGamma();
	double IRLS(int);
	double PIRLS(int);
    };

    class nlmerSp : public nlmer {
	rwSpFeMod fe;
    public:
	nlmerSp(Rcpp::S4 xp);
	void updateGamma();
	double IRLS(int);
	double PIRLS(int);
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
