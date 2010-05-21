// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_MER_H
#define LME4_MER_H

#include <Rcpp.h>
#include "MatrixNs.h"
#include "glmFamily.h"

namespace mer {
    class merResp;		// forward declaration

    class reModule {
	MatrixNs::chmFr d_L;
	MatrixNs::chmSp d_Lambda, d_Ut, d_Zt;
	Rcpp::IntegerVector Lind;
	Rcpp::NumericVector lower, theta;
	Rcpp::NumericVector d_u;
	double *d_ldL2, d_sqlLenU;
    public:
	reModule(Rcpp::S4);

	const Rcpp::NumericVector &u() const {return d_u;}
	const MatrixNs::chmFr &L() const {return d_L;}
	const MatrixNs::chmSp &Lambda() const {return d_Lambda;}
	const MatrixNs::chmSp &Ut() const {return d_Ut;}
	const MatrixNs::chmSp &Zt() const {return d_Zt;}
	double ldL2() const {return *d_ldL2;}
	double sqrLenU() const {return d_sqlLenU;}

	CHM_SP SupdateL(MatrixNs::chmSp const&) const;
	void DupdateL(MatrixNs::chmDn const&,MatrixNs::chmDn&) const;
	void rwUpdateL(Rcpp::NumericMatrix const&,
		       Rcpp::NumericVector const&,
		       Rcpp::NumericVector const&,
		       double*);
	void setU(std::vector<double> const&);
	void setU(Rcpp::NumericVector const&);
	void updateTheta(Rcpp::NumericVector const&);
	void updateU(merResp const&);
    };

    class merResp {
    protected:
	double *d_wrss;
	Rcpp::NumericVector d_offset, d_sqrtrwt, d_wtres, mu, weights, y;
    public:
	merResp(Rcpp::S4);
	void updateL(reModule const&);
	double updateMu(Rcpp::NumericVector const&);
	double updateWrss();
	double wrss() const {return *d_wrss;}
	const Rcpp::NumericVector &sqrtrwt() const {return d_sqrtrwt;}
	const Rcpp::NumericVector &wtres() const {return d_wtres;}
	const Rcpp::NumericVector &offset() const  {return d_offset;}

	Rcpp::NumericVector Utr, Vtr, cbeta, cu;
    };

    class feModule {
    protected:
	Rcpp::NumericVector d_beta;
	double *d_ldRX2;
    public:
	feModule(Rcpp::S4 xp);
	void setBeta(Rcpp::NumericVector const&);
	const Rcpp::NumericVector& beta() const;
	double ldRX2() const;
    };

    class deFeMod : public feModule {
    protected:
	MatrixNs::dgeMatrix d_X, d_RZX;
	MatrixNs::Cholesky d_RX;
    public:
	deFeMod(Rcpp::S4 xp);
	const MatrixNs::Cholesky &RX() const{return d_RX;}
	const MatrixNs::dgeMatrix &X() const{return d_X;}
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
	MatrixNs::chmSp d_X, d_RZX;
	MatrixNs::chmFr d_RX;
    public:
	spFeMod(Rcpp::S4 xp);
	const MatrixNs::chmSp& X() const {return d_X;}
	const MatrixNs::chmSp& RZX() const {return d_RZX;}
	const MatrixNs::chmFr& RX() const {return d_RX;}
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
	Rcpp::NumericMatrix d_sqrtXwt;
    public:
	rwResp(Rcpp::S4 xp);
	const Rcpp::NumericMatrix &sqrtXwt() const;
    };

    class glmerResp : public rwResp {
    protected:
	glmFamily family;
	Rcpp::NumericVector eta, muEta, n, d_var;
	void linkFun(){family.linkFun(eta, mu);}
	void linkInv(){family.linkInv(mu, eta);}
	void MuEta(){family.muEta(muEta, eta);}
	void variance(){family.variance(d_var, mu);}
    public:
	glmerResp(Rcpp::S4 xp);
	const Rcpp::NumericVector &var() const {return d_var;}
	Rcpp::NumericVector devResid(){return family.devResid(mu, weights, y);}
	double updateMu(Rcpp::NumericVector const &gamma);
	void updateSqrtRWt();
	void updateSqrtXWt();
    };
    
    class nlmerResp : public rwResp {
	Rcpp::Environment nlenv;
	Rcpp::Language nlmod;
	Rcpp::CharacterVector pnames;
    public:
	nlmerResp(Rcpp::S4 xp);
	double updateMu(Rcpp::NumericVector const &gamma);
    };
    
    class rwDeFeMod : public deFeMod {
	MatrixNs::dgeMatrix d_V;
    public:
	rwDeFeMod(Rcpp::S4 xp);
	void updateV(Rcpp::NumericMatrix const&);
	void updateRX(bool);
	void dpotrs(Rcpp::NumericVector&) const;
	const MatrixNs::dgeMatrix &V() const {return d_V;}
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
	double updateMu();
	double IRLS(int);
	double PIRLS(int);
    };

    class glmerSp : public glmer {
	rwSpFeMod fe;
    public:
	glmerSp(Rcpp::S4 xp);
	double updateMu();
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
	double updateMu();
	double IRLS(int);
	double PIRLS(int);
    };

    class nlmerSp : public nlmer {
	rwSpFeMod fe;
    public:
	nlmerSp(Rcpp::S4 xp);
	double updateMu();
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

	double pwrss() const;
	double deviance() const;
	double reCrit() const;
	double updateTheta(const Rcpp::NumericVector&);
    };
    
    template<typename T>
    inline lmer<T>::lmer(Rcpp::S4 xp) :
	re(Rcpp::S4(xp.slot("re"))),
	fe(Rcpp::S4(xp.slot("fe"))),
	resp(Rcpp::S4(xp.slot("resp")))
    {
	Rcpp::LogicalVector REML = xp.slot("REML");
	reml = (bool)*REML.begin();
    }

    template<typename T>
    double lmer<T>::l2PI = log(2. * PI);

/** 
 * Evaluate the penalized, weighted residual sum of squares
 * 
 * @return resp.wrss() + crossprod(re.u())
 */
    template<typename T> inline
    double lmer<T>::pwrss() const {
	const Rcpp::NumericVector u = re.u();
	return resp.wrss() +
	    std::inner_product(u.begin(), u.end(), u.begin(), double());
    }
/** 
 * Evaluate the deviance
 * 
 * @return ldL2 + n *(1 + log(2 * pi * pwrss/n))
 */
    template<typename T> inline
    double lmer<T>::deviance() const {
	double nn = (double)resp.wtres().size();
	return re.ldL2() + nn * (1 + l2PI + log(pwrss()/nn));
    }
/** 
 * Evaluate the REML criterion
 * 
 * @return ldL2 + ldRX2 + (n - p) * (1 + log(2 * pi * pwrss/(n - p)))
 */
    template<typename T> inline
    double lmer<T>::reCrit() const {
	double nmp = (double)(resp.wtres().size() - fe.beta().size());
	return re.ldL2()+fe.ldRX2()+nmp*(1 + l2PI + log(pwrss()/nmp));
    }

    template<typename T> inline
    double lmer<T>::updateTheta(Rcpp::NumericVector const &nt) {
	re.updateTheta(nt);
	resp.updateL(re);
	fe.updateRzxRx(re);
	fe.updateBeta(resp);
	re.updateU(resp);

	const Rcpp::NumericVector u = re.u(), offset = resp.offset();
	Rcpp::NumericVector b(u.size()), gamma(offset.size());
	std::copy(offset.begin(), offset.end(), gamma.begin());
	MatrixNs::chmDn bb(b), gg(gamma);

	re.Lambda().dmult('N', 1., 0., MatrixNs::chmDn(u), bb);
	re.Zt().dmult('T', 1., 1., bb, gg);
	fe.X().dmult('N', 1., 1., MatrixNs::chmDn(fe.beta()), gg);
	resp.updateMu(gamma);
	return reml ? reCrit() : deviance();
    }
}

#endif
