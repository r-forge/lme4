// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_RESPMODULE_H
#define LME4_RESPMODULE_H

#include <Rcpp.h>
#include "glmFamily.h"

namespace mer {
    class merResp {
    public:
	merResp(Rcpp::S4) throw (std::runtime_error);
	merResp(Rcpp::NumericVector y) throw (std::runtime_error);
	merResp(Rcpp::NumericVector y, Rcpp::NumericVector weights)
	    throw (std::runtime_error);
	merResp(Rcpp::NumericVector y, Rcpp::NumericVector weights,
	    Rcpp::NumericVector offset)
	    throw (std::runtime_error);

	void init();
	      Rcpp::NumericVector devResid() const;
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

    class lmerResp : public merResp {
    private:
	bool d_reml;
    public:
	lmerResp(Rcpp::S4) throw (std::runtime_error);
	lmerResp(bool,Rcpp::NumericVector) throw (std::runtime_error);
	lmerResp(bool,Rcpp::NumericVector,Rcpp::NumericVector)
	    throw (std::runtime_error);
	lmerResp(bool,Rcpp::NumericVector,Rcpp::NumericVector,
		 Rcpp::NumericVector) throw (std::runtime_error);

	const Rcpp::NumericVector&      mu() const {return d_mu;}
	const Rcpp::NumericVector&  offset() const {return d_offset;}
	const Rcpp::NumericMatrix& sqrtXwt() const {return d_sqrtXwt;}
	const Rcpp::NumericVector& sqrtrwt() const {return d_sqrtrwt;}
	const Rcpp::NumericVector& weights() const {return d_weights;}
	const Rcpp::NumericVector&   wtres() const {return d_wtres;}
	const Rcpp::NumericVector&       y() const {return d_y;}
	double                        wrss() const {return d_wrss;}
	double                   updateWts()       {return updateWrss();}
	double Laplace (double,double,double)const;
	double updateMu(const Rcpp::NumericVector&);
    };
	
    class glmerResp : public merResp {
    public:
	glmerResp(Rcpp::S4) throw (std::runtime_error);
	glmerResp(Rcpp::List, Rcpp::NumericVector y)
	    throw (std::runtime_error);
	glmerResp(Rcpp::List, Rcpp::NumericVector y,
		Rcpp::NumericVector wts) throw (std::runtime_error);
	glmerResp(Rcpp::List, Rcpp::NumericVector y,
		Rcpp::NumericVector wts,
		Rcpp::NumericVector offset) throw (std::runtime_error);
	glmerResp(Rcpp::List, Rcpp::NumericVector y,
		Rcpp::NumericVector wts,
		Rcpp::NumericVector offset,
		Rcpp::NumericVector n) throw (std::runtime_error);
	glmerResp(Rcpp::List, Rcpp::NumericVector y,
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
	double                       Laplace(double,double,double) const;
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
	glm::glmFamily d_fam;
	const Rcpp::NumericVector d_n;
	Rcpp::NumericVector d_eta;
    };

    class nlmerResp : public merResp {
	Rcpp::Environment d_nlenv;
	Rcpp::Language d_nlmod;
	Rcpp::CharacterVector d_pnames;
    public:
	nlmerResp(Rcpp::S4 xp);

	const Rcpp::NumericVector&        mu() const {return d_mu;}
	const Rcpp::NumericVector&    offset() const {return d_offset;}
	const Rcpp::NumericMatrix&   sqrtXwt() const {return d_sqrtXwt;}
	const Rcpp::NumericVector&   sqrtrwt() const {return d_sqrtrwt;}
	const Rcpp::NumericVector&   weights() const {return d_weights;}
	const Rcpp::NumericVector&     wtres() const {return d_wtres;}
	const Rcpp::NumericVector&         y() const {return d_y;}
	double Laplace (double, double, double) const;
	double updateMu(const Rcpp::NumericVector&) throw(std::runtime_error);
    };
    
}

#endif
