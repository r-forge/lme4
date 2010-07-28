// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
#ifndef LME4_UTILITIES_H
#define LME4_UTILITIES_H

#include <Rcpp.h>
#include "MatrixNs.h"

namespace mer {
    double compareVecWt(const Rcpp::NumericVector&,
			const Rcpp::NumericVector&,
			const Rcpp::NumericVector&);
    void showCHM_DN    (const_CHM_DN, const std::string&);
    void showCHM_FR    (const_CHM_FR, const std::string&);
    void showCHM_SP    (const_CHM_SP, const std::string&);
    void showCHM_SP    (const MatrixNs::chmSp&, const std::string&);
    void showdMat      (const MatrixNs::ddenseMatrix&, const std::string&);
    void showdbl       (const double*, const char*, int);
    void showdbl       (Rcpp::NumericVector const&, const char*);
    void showincr      (double, double, double,
		        const Rcpp::NumericVector&, const char*);
    void showint       (const int*, const char*, int);

    // Add the fixed effects parameters and the random effects values
    // as attributes of the deviance or REML criterion.  The purpose
    // is to be able to treat the deviance or REML as the objective in
    // an optimization but also retain some of the status information.
    Rcpp::NumericVector mkans(double,
			      const Rcpp::NumericVector&,
			      const Rcpp::NumericVector&,
			      double,double,double);

}

#endif 
