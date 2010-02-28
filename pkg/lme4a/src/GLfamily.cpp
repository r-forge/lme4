#include "GLfamily.hpp"
#include "lme4utils.hpp"

static clogloglink clogloglnk;
static identitylink identitylnk;
static logitlink logitlnk;
static loglink loglnk;
static probitlink probitlnk;
static sqrtlink sqrtlnk;
static constvar constvr;
static mu1muvar mu1muvr;
static mu2var mu2vr;
static mu3var mu3vr;
static muvar muvr;

void GLfamily::initGL(SEXP rho) {
    static std::map<std::string, GLlink*> lpts;
    static std::map<std::string, GLvar*> vpts;
    if (!lpts.count("identity")) { // initialize the map contents
	lpts["identity"] = &identitylnk;
	lpts["logit"] = &logitlnk;
	lpts["probit"] = &probitlnk;
	lpts["sqrt"] = &probitlnk;
	vpts["gaussian"] = &constvr;
	vpts["binomial"] = &mu1muvr;
	vpts["Gamma"] = &mu2vr;
	vpts["inverse.gaussian"] = &mu3vr;
	vpts["poisson"] = &muvr;
    }
    
    family = findVarBound(rho, lme4_familySym);
    SEXP nms = getAttrib(family, R_NamesSymbol);
    if (!isNewList(family) ||!isString(nms) ||
	LENGTH(nms) != LENGTH(family))
	error(_("Object \"%s\" must be a named list"), "family");
    fname = CHAR(STRING_ELT(getListElement(family, nms, "family"), 0));
    lname = CHAR(STRING_ELT(getListElement(family, nms, "link"), 0));
				// assign the link and variance class
    if (lpts.count(lname)) lnk = lpts[lname];
    else error(_("no compiled link function named \"%s\""), lname);
    if (vpts.count(fname)) var = vpts[fname];
    else error(_("no compiled variance function for \"%s\" distribution"),
	       fname);
}    

const double GLlink::LTHRESH = 30;
const double GLlink::MLTHRESH = -30;
const double GLlink::MPTHRESH = qnorm5(DOUBLE_EPS, 0, 1, 1, 0);
const double GLlink::PTHRESH = -MPTHRESH;
const double GLlink::INVEPS = 1/DOUBLE_EPS;

void logitlink::link(double *eta, const double* mu, int n) {
    for (int i = 0; i < n; i++) {
	double mui = mu[i];
	eta[i] = log(mui/(1 - mui));
    }
}

void logitlink::linkinv(double *mu, double *muEta,
			const double* eta, int n) {
    for (int i = 0; i < n; i++) {
	double etai = eta[i], tmp;
	tmp = (etai < MLTHRESH) ? DOUBLE_EPS : ((etai > LTHRESH) ?
						INVEPS : exp(etai));
	mu[i] = tmp/(1 + tmp);
	muEta[i] = mu[i] * (1 - mu[i]);
    }
}

void probitlink::link(double *eta, const double* mu, int n) {
    for (int i = 0; i < n; i++) eta[i] = qnorm5(mu[i], 0, 1, 1, 0);
}

void probitlink::linkinv(double *mu, double *muEta,
			 const double* eta, int n)
{
    for (int i = 0; i < n; i++) {
	double etai = eta[i], tmp;
	mu[i] = (etai < MPTHRESH) ? DOUBLE_EPS :
	    ((etai > PTHRESH) ? 1 - DOUBLE_EPS : pnorm5(etai, 0, 1, 1, 0));
	tmp = dnorm4(etai, 0, 1, 0);
	muEta[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
    }
}

void clogloglink::link(double *eta, const double* mu, int n) {
    for (int i = 0; i < n; i++) eta[i] = log(-log(1 - mu[i]));
}

void clogloglink::linkinv(double *mu, double *muEta,
			  const double* eta, int n)
{
    for (int i = 0; i < n; i++) {
	double etai = eta[i], tmp, t2;
	tmp = (etai < MLTHRESH) ? DOUBLE_EPS : ((etai > LTHRESH) ?
						INVEPS : exp(etai));
	t2 = -expm1(-tmp);
	mu[i] = (t2 < DOUBLE_EPS) ? DOUBLE_EPS :
	    (t2 > 1 - DOUBLE_EPS ? 1 - DOUBLE_EPS : t2);
	muEta[i] = tmp * exp(-tmp);
    }
}

void identitylink::link(double *eta, const double* mu, int n) {
    for (int i = 0; i < n; i++) eta[i] = mu[i];
}

void identitylink::linkinv(double *mu, double *muEta,
			  const double* eta, int n)
{
    for (int i = 0; i < n; i++) {
	mu[i] = eta[i];
	muEta[i] = 1;
    }
}

void loglink::link(double *eta, const double* mu, int n) {
    for (int i = 0; i < n; i++) eta[i] = log(mu[i]);
}

void loglink::linkinv(double *mu, double *muEta, const double* eta, int n)
{
    for (int i = 0; i < n; i++) {
	double tmp = exp(eta[i]);
	muEta[i] = mu[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
    }
}

void sqrtlink::link(double *eta, const double* mu, int n) {
    for (int i = 0; i < n; i++) eta[i] = sqrt(mu[i]);
}

void sqrtlink::linkinv(double *mu, double *muEta, const double* eta, int n)
{
    for (int i = 0; i < n; i++) {
	double etai = eta[i];
	mu[i] = etai * etai;
	muEta[i] = 2 * etai;
    }
}

// void GLvar::devResid(double *ans, const double *mu,
// 			const double *pw, const double *y,
// 			const int *fac, int n) {
//     error(_("Should not be called"));
// }

// void GLvar::varFunc(double *var, const double *mu, int n) {
//     error(_("Should not be called"));
// }


void constvar::devResid(double *ans, const double *mu,
			const double *pw, const double *y,
			const int *fac, int n) {
    for (int i = 0; i < n; i++) {
	double ri = y[i] - mu[i];
	ans[fac ? (fac[i] - 1) : 0] += (pw ? pw[i] : 1) * ri * ri;
    }
}

void muvar::devResid(double *ans, const double *mu,
		     const double *pw, const double *y,
		     const int *fac, int n)
{
    for (int i = 0; i < n; i++) {
	ans[fac ? (fac[i] - 1) : 0] += 2 * (pw ? pw[i] : 1) *
	    (y_log_y(y[i], mu[i]) - (y[i] - mu[i]));
    }
}

void mu1muvar::devResid(double *ans, const double *mu,
			const double *pw, const double *y,
			const int *fac, int n)
{
    for (int i = 0; i < n; i++) {
	ans[fac ? (fac[i] - 1) : 0] += 2 * (pw ? pw[i] : 1) *
	    (y_log_y(y[i], mu[i]) + y_log_y(1 - y[i], 1 - mu[i]));

    }
}

void mu2var::devResid(double *ans, const double *mu,
			const double *pw, const double *y,
			const int *fac, int n)
{
    for (int i = 0; i < n; i++) {
	ans[fac ? (fac[i] - 1) : 0] += 2 * (pw ? pw[i] : 1) *
	    (y_log_y(y[i], mu[i]) - (y[i] - mu[i])/mu[i]);
    }
}

void mu3var::devResid(double *ans, const double *mu,
		      const double *pw, const double *y,
		      const int *fac, int n)
{
    for (int i = 0; i < n; i++) {
	double ri = y[i] - mu[i];
	ans[fac ? (fac[i] - 1) : 0] += 2 * (pw ? pw[i] : 1) *
	    (ri * ri)/(y[i] * mu[i] * mu[i]);
    }
}

void constvar::varFunc(double *var, const double *mu, int n)
{
    for (int i = 0; i < n; i++)
	var[i] = 1.;
}

void mu1muvar::varFunc(double *var, const double *mu, int n)
{
    for (int i = 0; i < n; i++) {
	double mui = mu[i];
	if (mui <= 0 || mui >= 1)
	    error(_("mu[%d] = %g is not in (0,1)"), i + 1, mui);
	var[i] = mui * (1 - mui);
    }
}

void muvar::varFunc(double *var, const double *mu, int n)
{
    for (int i = 0; i < n; i++) {
	double mui = mu[i];
	if (mui <= 0)
	    error(_("mu[%d] = %g is not positive"), i + 1, mui);
	var[i] = mui;
    }
}

void mu2var::varFunc(double *var, const double *mu, int n)
{
    for (int i = 0; i < n; i++) {
	double mui = mu[i];
	if (mui <= 0)
	    error(_("mu[%d] = %g is not positive"), i + 1, mui);
	var[i] = mui * mui;
    }
}

void mu3var::varFunc(double *var, const double *mu, int n)
{
    for (int i = 0; i < n; i++) {
	double mui = mu[i];
	if (mui <= 0)
	    error(_("mu[%d] = %g is not positive"), i + 1, mui);
	var[i] = mui * mui * mui;
    }
}

