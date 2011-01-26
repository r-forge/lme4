#include "merMod.h"
#include <R_ext/BLAS.h>

using namespace Rcpp;
using namespace MatrixNs;
using namespace std;
    
/** 
 * Evaluate the profiled deviance or the profiled REML criterion from
 * a theta parameter value for an lmer model.
 * 
 * @param xp merMod object
 * @param theta new value of variance component parameters
 * @param beta value of fixed-effects parameters (ignored unless alg == 2)
 * @param u0 starting value for random effects (ignored unless alg == 2 or 3)
 * @param verb level of verbosity, verb == 2 is very verbose
 * @param alg alg == 1 => IRLS, alg == 2 => PIRLS, alg == 3 => PIRLSBeta
 * 
 * @return a deviance evaluation
 */

extern "C" 
SEXP merDeviance(SEXP xpp, SEXP thetap, SEXP betap, SEXP u0p, SEXP verbp, SEXP algp) {
    try {
	S4 xp(xpp);
	NumericVector theta(thetap), beta(betap), u0(u0p);
	int verb = as<int>(verbp), alg = as<int>(algp);
	S4 re(xp.slot("re")), fe(xp.slot("fe")), resp(xp.slot("resp"));
	bool de = fe.is("deFeMod");
	if (!de && !fe.is("spFeMod"))
	    throw runtime_error("fe slot is neither deFeMod nor spFeMod");
	int rtype = resp.is("lmerResp") ? 1 : (resp.is("glmRespMod") ? 2 :
					       (resp.is("nlsRespMod") ? 3 : 0));
	if (rtype == 0) throw runtime_error("unknown respMod type in merDeviance");
	
	if (alg < 1 || alg > 3) throw range_error("alg must be 1, 2 or 3");
	mer::Alg aa = (alg == 1) ? mer::Beta : ((alg == 2) ? mer::U : mer::BetaU);
	
	switch(rtype) {
	case 1:
	    if (de) {
		mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
		return lm.LMMdeviance(theta);
	    } else {
		mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
		return lm.LMMdeviance(theta);
	    }
	    break;
	case 2:
	    if (de) {
		mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
		return glmr.PIRLS(theta, beta, u0, verb, aa);
	    } else {
		mer::mer<mer::spFeMod,mer::glmerResp> glmr(xp);
		return glmr.PIRLS(theta, beta, u0, verb, aa);
	    }
	    break;
	case 3:
	    if (de) {
		mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
		return nlmr.PIRLS(theta, beta, u0, verb, aa);
	    } else {
		mer::mer<mer::spFeMod,mer::nlmerResp> nlmr(xp);
		return nlmr.PIRLS(theta, beta, u0, verb, aa);
	    }
	}
    } catch( std::exception& __ex__ ) {
	forward_exception_to_r( __ex__ );
    } catch(...) {
	::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}

extern "C"
SEXP updateDc(SEXP xpp, SEXP thp, SEXP betap, SEXP up) {
    try {
	S4 xp(xpp);
	NumericVector th(thp), beta(betap), u(up);
	S4 fe(xp.slot("fe")), resp(xp.slot("resp"));
	bool de(fe.is("deFeMod"));
	if (!de && !fe.is("spFeMod"))
	    throw runtime_error("fe slot is neither deFeMod nor spFeMod");
	if (resp.is("lmerResp")) {
	    if (de) {
		mer::mer<mer::deFeMod,mer::lmerResp> lm(xp);
		return lm.updateDcmp(th, beta, u);
	    } else {
		mer::mer<mer::spFeMod,mer::lmerResp> lm(xp);
		return lm.updateDcmp(th, beta, u);
	    }
	}
	if (resp.is("glmRespMod")) {
	    if (de) {
		mer::mer<mer::deFeMod,mer::glmerResp> glmr(xp);
		return glmr.updateDcmp(th, beta, u);
	    } else {
		mer::mer<mer::spFeMod,mer::glmerResp> glmr(xp);
		return glmr.updateDcmp(th, beta, u);
	    }
	} else if (resp.is("nlsRespMod")) {
	    if (de) {
		mer::mer<mer::deFeMod,mer::nlmerResp> nlmr(xp);
		return nlmr.updateDcmp(th, beta, u);
	    } else {
		mer::mer<mer::spFeMod,mer::nlmerResp> nlmr(xp);
		return nlmr.updateDcmp(th, beta, u);
	    }
	}
	throw runtime_error("resp slot is not lmerResp or glmRespMod or nlsRespMod");
    } catch( std::exception& __ex__ ) {
	forward_exception_to_r( __ex__ );
    } catch(...) {
	::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}    

extern "C" 
SEXP reTrmsCondVar(SEXP xpp, SEXP scalep) {
    try {
	S4 xp(xpp);
	double scale = as<double>(scalep);
	mer::reTrms trms(xp);
	return trms.condVar(scale);
    } catch( std::exception& __ex__ ) {
	forward_exception_to_r( __ex__ );
    } catch(...) {
	::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}
    
extern "C" 
SEXP glmIRLS(SEXP xpp, SEXP verbp) {
    try {
	S4 xp(xpp);
	int verb = as<int>(verbp);
	S4 pm = xp.slot("pred");
	if (pm.is("dPredModule")) {
	    glm::mod<matMod::dPredModule,mer::glmerResp> m(xp);
	    return m.IRLS(verb);
	} else if(pm.is("sPredModule")) {
	    glm::mod<matMod::sPredModule,mer::glmerResp> m(xp);
	    return m.IRLS(verb);
	} else throw runtime_error("Unknown linear predictor module type");
    } catch( std::exception& __ex__ ) {
	forward_exception_to_r( __ex__ );
    } catch(...) {
	::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}
