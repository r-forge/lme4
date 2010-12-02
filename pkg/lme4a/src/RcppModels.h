// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "glmFamily.h"

namespace glm {
    class lmResp {
    public:
	lmResp() {}
	lmResp(Rcpp::NumericVector y) throw (std::runtime_error);
	lmResp(Rcpp::NumericVector y, Rcpp::NumericVector weights)
	    throw (std::runtime_error);
	lmResp(Rcpp::NumericVector y, Rcpp::NumericVector weights,
	    Rcpp::NumericVector offset)
	    throw (std::runtime_error);

	void init();
	const Rcpp::NumericVector&      mu() const {return d_mu;}
	const Rcpp::NumericVector&  offset() const {return d_offset;}
	const Rcpp::NumericMatrix& sqrtXwt() const {return d_sqrtXwt;}
	const Rcpp::NumericVector& sqrtrwt() const {return d_sqrtrwt;}
	const Rcpp::NumericVector& weights() const {return d_weights;}
	const Rcpp::NumericVector&   wtres() const {return d_wtres;}
	const Rcpp::NumericVector&       y() const {return d_y;}
	double                        wrss() const {return d_wrss;}
	double                    updateMu(const Rcpp::NumericVector&);
	double                   updateWts()       {return updateWrss();}
	double                  updateWrss();
    protected:
	double                     d_wrss;
	const Rcpp::NumericVector  d_y, d_weights, d_offset;
	Rcpp::NumericVector        d_mu, d_sqrtrwt, d_wtres;
	Rcpp::NumericMatrix        d_sqrtXwt;
    };
	
    class glmResp : public lmResp {
    public:
	glmResp(Rcpp::List, Rcpp::NumericVector y)
	    throw (std::runtime_error);
	glmResp(Rcpp::List, Rcpp::NumericVector y,
		Rcpp::NumericVector wts) throw (std::runtime_error);
	glmResp(Rcpp::List, Rcpp::NumericVector y,
		Rcpp::NumericVector wts,
		Rcpp::NumericVector offset) throw (std::runtime_error);
	glmResp(Rcpp::List, Rcpp::NumericVector y,
		Rcpp::NumericVector wts,
		Rcpp::NumericVector offset,
		Rcpp::NumericVector n) throw (std::runtime_error);
	glmResp(Rcpp::List, Rcpp::NumericVector y,
		Rcpp::NumericVector wts,
		Rcpp::NumericVector offset,
		Rcpp::NumericVector n,
		Rcpp::NumericVector eta) throw (std::runtime_error);

	const Rcpp::NumericVector&       eta() const {return d_eta;}
	const Rcpp::NumericVector&        mu() const {return d_mu;}
	const Rcpp::NumericVector&    offset() const {return d_offset;}
	const Rcpp::NumericMatrix&   sqrtXwt() const {return d_sqrtXwt;}
	const Rcpp::NumericVector&   sqrtrwt() const {return d_sqrtrwt;}
	const Rcpp::NumericVector&   weights() const {return d_weights;}
	const Rcpp::NumericVector&     wtres() const {return d_wtres;}
	const Rcpp::NumericVector&         y() const {return d_y;}
	const std::string&            family() const {return d_fam.fam();}
	const std::string&              link() const {return d_fam.lnk();}
//	double                      deviance() const;
	double                 residDeviance() const;
	double                      updateMu(const Rcpp::NumericVector&);
	double                     updateWts();
	double                          wrss() const {return d_wrss;}
	Rcpp::NumericVector         devResid() const;
	Rcpp::NumericVector            muEta() const;
        Rcpp::NumericVector        sqrtWrkWt() const;
	Rcpp::NumericVector         variance() const;
        Rcpp::NumericVector        wrkResids() const;
        Rcpp::NumericVector          wrkResp() const;
    protected:
	glmFamily d_fam;
	const Rcpp::NumericVector d_n;
	Rcpp::NumericVector d_eta;
    };

    class predMod {
    protected:
	Rcpp::NumericVector d_coef, d_coef0, d_delta;
	arma::vec a_coef, a_coef0, a_delta;
    public:
	predMod(int);
	predMod(Rcpp::NumericVector);

	const Rcpp::NumericVector& delta() const {return d_delta;}
	const Rcpp::NumericVector& coef0() const {return d_coef0;}

	void installCoef0();
	void setCoef0(const Rcpp::NumericVector&)
	    throw (std::runtime_error);
    };

    class densePred : predMod {
    protected:
	Rcpp::NumericMatrix d_X;
	arma::mat a_X;
    public:
	densePred(Rcpp::NumericMatrix) throw (std::runtime_error);
	densePred(Rcpp::NumericMatrix, Rcpp::NumericVector)
	    throw (std::runtime_error);

	const Rcpp::NumericMatrix&     X() const {return d_X;}
	const arma::mat&              aX() const {return a_X;}	
	const Rcpp::NumericVector&  coef() const {return d_coef;}
	const Rcpp::NumericVector& coef0() const {return d_coef0;}
	const Rcpp::NumericVector& delta() const {return d_delta;}
	void installCoef0();
	void setCoef0(const Rcpp::NumericVector&)
	    throw (std::runtime_error);
// Choose better names than gamma and dgamma.  These conflict with
// gamma functions, gamma distributions, the Gamma family, ...
	Rcpp::NumericVector gamma(const Rcpp::NumericVector&,
				  const Rcpp::NumericVector&)
	    throw (std::runtime_error);
	Rcpp::NumericVector dgamma(double fac);
    };    
}
