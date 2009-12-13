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

//#define USING_MAP_TEMPLATE
#ifdef USING_MAP_TEMPLATE
#include <string>
#include <map>
#else
static const int nlink = 6;
static const int nvar = 5;
#endif

void GLfamily::initGL(SEXP rho) {
#ifdef USING_MAP_TEMPLATE
    std::map<string, GLlink> links;
    std::map<string, GLvar> vars;
    links["cloglog"] = clogloglnk;
    links["identity"] = identitylnk;
    links["logit"] = logitlnk;
    links["probit"] = probitlnk;
    links["sqrt"] = probitlnk;
    vars["gaussian"] = constvr;
    vars["binomial"] = mu1muvr;
    vars["Gamma"] = mu2vr;
    vars["poisson"] = muvr;
    vars["inverse.gaussian"] = mu3vr;
#else    
    GLlink **links = new GLlink*[nlink];
    links[0] = &clogloglnk;
    links[1] = &identitylnk;
    links[2] = &logitlnk;
    links[3] = &loglnk;
    links[4] = &probitlnk;
    links[5] = &sqrtlnk;
    GLvar **vars = new GLvar*[nvar];
    vars[0] = &constvr;
    vars[1] = &mu1muvr;
    vars[2] = &mu2vr;
    vars[3] = &mu3vr;
    vars[4] = &muvr;
#endif
    
    family = findVarBound(rho, install("family"));
    SEXP nms = getAttrib(family, R_NamesSymbol);
    if (!isNewList(family) ||!isString(nms) ||
	LENGTH(nms) != LENGTH(family))
	error(_("Object \"%s\" must be a named list"), "family");
    fname = CHAR(STRING_ELT(getListElement(family, nms, "family"), 0));
    lname = CHAR(STRING_ELT(getListElement(family, nms, "link"), 0));
    Rprintf("family name: %s\n", fname);
    Rprintf("link name in family: %s\n", lname);
    // assign the link and variance function by comparing names
#ifdef USING_MAP_TEMPLATE
    lnk = links[lname];
    var = vars[fname];
#else
    for (int i = 0; i < nlink; i++)
	if (!strcmp(lname, links[i]->nm())) lnk = links[i];
    for (int i = 0; i < nvar; i++)
	if (!strcmp(fname, vars[i]->distnm())) var = vars[i];
#endif
    Rprintf("variance name: %s\n", var->distnm());
    Rprintf("link name: %s\n", lnk->nm());
#ifndef USING_MAP_TEMPLATE
    delete[] links;
    delete[] vars;
#endif
}    

const double GLlink::LTHRESH = 30;
const double GLlink::MLTHRESH = -30;
const double GLlink::MPTHRESH = qnorm5(DOUBLE_EPS, 0, 1, 1, 0);
const double GLlink::PTHRESH = -MPTHRESH;
const double GLlink::INVEPS = 1/DOUBLE_EPS;

void logitlink::linkinv(double *mu, double *muEta,
			const double* eta, int n)
{
    for (int i = 0; i < n; i++) {
	double etai = eta[i], tmp;
	tmp = (etai < MLTHRESH) ? DOUBLE_EPS : ((etai > LTHRESH) ?
						INVEPS : exp(etai));
	mu[i] = tmp/(1 + tmp);
	muEta[i] = mu[i] * (1 - mu[i]);
    }
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

void identitylink::linkinv(double *mu, double *muEta,
			  const double* eta, int n)
{
    for (int i = 0; i < n; i++) {
	mu[i] = eta[i];
	muEta[i] = 1;
    }
}

void loglink::linkinv(double *mu, double *muEta, const double* eta, int n)
{
    for (int i = 0; i < n; i++) {
	double tmp = exp(eta[i]);
	muEta[i] = mu[i] = (tmp < DOUBLE_EPS) ? DOUBLE_EPS : tmp;
    }
}

void sqrtlink::linkinv(double *mu, double *muEta, const double* eta, int n)
{
    for (int i = 0; i < n; i++) {
	double etai = eta[i];
	mu[i] = etai * etai;
	muEta[i] = 2 * etai;
    }
}

void constvar::devResid(double *ans, const double *mu,
			const double *pw, const double *y,
			const int *fac, int n)
{
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

