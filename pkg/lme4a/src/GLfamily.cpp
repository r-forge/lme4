#include "GLfamily.hpp"

#define NO_C_HEADERS
#include <cstdio>
#include <cmath>
#include <cstring>
extern "C" {
#include <R.h>
#include <Rmath.h>
}

#ifdef ENABLE_NLS	    // Allow for translation of error messages
#include <libintl.h>
#define _(String) dgettext ("lme4", String)
#else
#define _(String) (String)
#endif

const double GLlink::LTHRESH = 30;
const double GLlink::MLTHRESH = -30;
const double GLlink::MPTHRESH = qnorm5(DOUBLE_EPS, 0, 1, 1, 0);
const double GLlink::PTHRESH = -MPTHRESH;
const double GLlink::INVEPS = 1/DOUBLE_EPS;

/**
 * Evaluate y * log(y/mu) with the correct limiting value at y = 0.
 *
 * @param y 
 * @param mu
 *
 * @return y * log(y/mu) for y > 0, 0 for y == 0.
 */
inline double y_log_y(double y, double mu)
{
    return (y) ? (y * log(y/mu)) : 0;
}

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

