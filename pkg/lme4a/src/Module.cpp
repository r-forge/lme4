#include "glmFamily.h"
#include "RcppModels.h"

using namespace Rcpp;

RCPP_MODULE(lme4) {

using namespace glm;
    
class_<glmFamily>("glmFamily")

    .constructor<List>()

    .property("family", &glmFamily::fam)
    .method("linkFun",  &glmFamily::linkFun)
    .method("linkInv",  &glmFamily::linkInv)
    .method("muEta",    &glmFamily::muEta)
    .method("devResid", &glmFamily::devResid)
    .method("variance", &glmFamily::variance)
    ;

class_<lmResp>("lmResp")

    .constructor<NumericVector>()
    .constructor<NumericVector,NumericVector>()
    .constructor<NumericVector,NumericVector,NumericVector>()

    .property("mu",      &lmResp::mu)
    .property("offset",  &lmResp::offset)
    .property("sqrtXwt", &lmResp::sqrtXwt)
    .property("sqrtrwt", &lmResp::sqrtrwt)
    .property("weights", &lmResp::weights)
    .property("wrss",    &lmResp::wrss)
    .property("wtres",   &lmResp::wtres)
    .property("y",       &lmResp::y)

    .method("updateMu",  &lmResp::updateMu)
    .method("updateWts", &lmResp::updateWts)
    ;

class_<glmResp>("glmResp")

    .constructor<List,NumericVector>()
    .constructor<List,NumericVector,NumericVector>()
    .constructor<List,NumericVector,NumericVector,NumericVector>()
    .constructor<List,NumericVector,NumericVector,NumericVector,NumericVector>()
 
    .property("devResid",    &glmResp::devResid,
	    "deviance residuals")
    .property("eta",         &glmResp::eta,
	      "current value of the linear predictor")
    .property("family",      &glmResp::family,
	      "name of the glm family")
    .property("link",        &glmResp::link,
	      "name of the glm link")
    .property("mu",          &glmResp::mu,
	      "current value of the mean vector")
    .property("muEta",       &glmResp::muEta,
	      "current value of d mu/d eta")
    .property("offset",      &glmResp::offset,
	      "offset vector - const")
    .property("residDeviance", &glmResp::residDeviance,
	      "sum of the deviance residuals")
    .property("sqrtXwt",     &glmResp::sqrtXwt,
	      "square root of the weights applied to the model matrix")
    .property("sqrtrwt",     &glmResp::sqrtrwt,
	      "square root of the weights applied to the residuals")
    .property("sqrtWrkWt",     &glmResp::sqrtrwt,
	      "square root of the weights applied to the working response")
    .property("weights",     &glmResp::weights,
	      "prior weights - const")
    .property("variance",    &glmResp::variance,
	      "current (unscaled) variances")
    .property("wtres",       &glmResp::wtres,
	      "current value of the weighted residuals")
    .property("wrss",        &glmResp::wrss,
	      "weighted residual sum of squares")
    .property("wrkResids",   &glmResp::wrkResids,
	      "working residuals - on the eta scale")
    .property("wrkResp",     &glmResp::wrkResp,
	      "working response - on the eta scale")
    .property("y",           &glmResp::y,
	      "numeric response vector - const")

    .method("updateMu",      &glmResp::updateMu,
	    "update mu and derived quantities from a new value of eta")
    .method("updateWts",     &glmResp::updateWts,
	    "update the residual and X weights.")
    ;

class_<densePred>("densePred")

    .constructor<NumericMatrix>()
    .constructor<NumericMatrix, NumericVector>()

    .property("coef",  &densePred::coef,
	      "last coefficient vector evaluated")
    .property("coef0", &densePred::coef0, &densePred::setCoef0,
	      "base coefficient vector")
    .property("delta", &densePred::delta,
	      "increment in coefficient vector")
    .property("X",     &densePred::X,
	      "dense model matrix")

    .method("gamma",   &densePred::gamma,
	    "evaluate delta from Xwts and wtres; return eta for full step factor")
    .method("dgamma",  &densePred::dgamma,
	    "return eta for specified step factor")
    .method("installCoef0", &densePred::installCoef0,
	    "install the current coef as coef0 for the next step")
    ;

}
