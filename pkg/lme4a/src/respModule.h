// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_RESPMODULE_H
#define LME4_RESPMODULE_H

#include <Rcpp.h>
#include "glmFamily.h"

namespace mer {
    class merResp {
    protected:
	Rcpp::S4             d_xp;
	double               d_wrss;
	Rcpp::NumericVector  d_offset, d_weights, d_y; // read-only slots
	Rcpp::NumericVector  d_sqrtrwt, d_wtres, d_mu; // writeable slots
	Rcpp::NumericMatrix  d_sqrtXwt;                // writeable
    public:
	merResp(Rcpp::S4);

	Rcpp::NumericVector       devResid() const;

	const Rcpp::NumericVector&      mu() const {return d_mu;}
	const Rcpp::NumericVector&  offset() const {return d_offset;}
	const Rcpp::NumericMatrix& sqrtXwt() const {return d_sqrtXwt;}
	const Rcpp::NumericVector& sqrtrwt() const {return d_sqrtrwt;}
	const Rcpp::NumericVector&   wtres() const {return d_wtres;}
	const Rcpp::NumericVector&       y() const {return d_y;}
	const Rcpp::S4&                 xp() const {return d_xp;}
	double                        wrss() const {return d_wrss;}
	double                   updateWts()       {return updateWrss();}
	double                  updateWrss();
    };

    class lmerResp : public merResp {
	int d_reml;
    public:
	lmerResp(Rcpp::S4);

	double Laplace (double,double,double)const;
	double updateMu(const Rcpp::NumericVector&);
    };

    class glmerResp : public merResp {
    protected:
	glm::glmFamily           family;
	Rcpp::NumericVector  d_eta, d_n;
    public:
	glmerResp(Rcpp::S4);

	Rcpp::NumericVector   devResid() const;

	const Rcpp::NumericVector& eta() const {return d_eta;}

	double                 Laplace(double,double,double) const;
	double                updateMu(const Rcpp::NumericVector&);
	double               updateWts();
    };
    
    class nlmerResp : public merResp {
	Rcpp::Environment nlenv;
	Rcpp::Language nlmod;
	Rcpp::CharacterVector pnames;
    public:
	nlmerResp(Rcpp::S4 xp);

	double Laplace (double, double, double) const;
	double updateMu(const Rcpp::NumericVector&);
    };
    
}

#endif
