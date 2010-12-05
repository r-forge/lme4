#include "glmFamily.h"
#include "respModule.h"

using namespace Rcpp;

RCPP_MODULE(lme4) {

class_<glm::glmFamily>("glmFamily")

    .constructor<List>()

    .property("family", &glm::glmFamily::fam)
    .property("link",   &glm::glmFamily::lnk)
    .method("linkFun",  &glm::glmFamily::linkFun)
    .method("linkInv",  &glm::glmFamily::linkInv)
    .method("muEta",    &glm::glmFamily::muEta)
    .method("devResid", &glm::glmFamily::devResid)
    .method("variance", &glm::glmFamily::variance)
    ;

using namespace mer;

class_<lmerResp>("lmerResp")

    .constructor<S4>()
    .constructor<int,NumericVector>()
    .constructor<int,NumericVector,NumericVector>()
    .constructor<int,NumericVector,NumericVector,NumericVector>()

//    .property("devResid", &lmerResp::devResid)
    .property("mu",       &lmerResp::mu)
    .property("offset",   &lmerResp::offset)
    .property("sqrtXwt",  &lmerResp::sqrtXwt)
    .property("sqrtrwt",  &lmerResp::sqrtrwt)
    .property("weights",  &lmerResp::weights)
    .property("wrss",     &lmerResp::wrss)
    .property("wtres",    &lmerResp::wtres)
    .property("y",        &lmerResp::y)

    .method("updateMu",   &lmerResp::updateMu)
    .method("updateWts",  &lmerResp::updateWts)
    ;

class_<glmerResp>("glmerResp")

    .constructor<S4>()
    .constructor<List,NumericVector>()
    .constructor<List,NumericVector,NumericVector>()
    .constructor<List,NumericVector,NumericVector,NumericVector>()
    .constructor<List,NumericVector,NumericVector,NumericVector,NumericVector>()
 
    .property("devResid",      &glmerResp::devResid,
	    "deviance residuals")
    .property("eta",           &glmerResp::eta,
	      "current value of the linear predictor")
    .property("family",        &glmerResp::family,
	      "name of the glm family")
    .property("link",          &glmerResp::link,
	      "name of the glm link")
    .property("mu",            &glmerResp::mu,
	      "current value of the mean vector")
    .property("muEta",         &glmerResp::muEta,
	      "current value of d mu/d eta")
    .property("offset",        &glmerResp::offset,
	      "offset vector - const")
    .property("residDeviance", &glmerResp::residDeviance,
	      "sum of the deviance residuals")
    .property("sqrtXwt",       &glmerResp::sqrtXwt,
	      "square root of the weights applied to the model matrix")
    .property("sqrtrwt",       &glmerResp::sqrtrwt,
	      "square root of the weights applied to the residuals")
    .property("sqrtWrkWt",     &glmerResp::sqrtrwt,
	      "square root of the weights applied to the working response")
    .property("weights",       &glmerResp::weights,
	      "prior weights - const")
    .property("variance",      &glmerResp::variance,
	      "current (unscaled) variances")
    .property("wtres",         &glmerResp::wtres,
	      "current value of the weighted residuals")
    .property("wrss",          &glmerResp::wrss,
	      "weighted residual sum of squares")
    .property("wrkResids",     &glmerResp::wrkResids,
	      "working residuals - on the eta scale")
    .property("wrkResp",       &glmerResp::wrkResp,
	      "working response - on the eta scale")
    .property("y",             &glmerResp::y,
	      "numeric response vector - const")

    .method("updateMu",      &glmerResp::updateMu,
	    "update mu and derived quantities from a new value of eta")
    .method("updateWts",     &glmerResp::updateWts,
	    "update the residual and X weights.")
    ;
}
