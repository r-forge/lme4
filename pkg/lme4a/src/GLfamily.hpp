#ifndef LME4_GLFAMILY_HPP
#define LME4_GLFAMILY_HPP

#include "lme4utils.h"		// for definition of SEXP

class GLlink {
public:
    virtual void linkinv(double *mu, double *muEta,
			 const double* eta, int n) = 0;
    virtual const char* nm() = 0;
    static const double LTHRESH, MLTHRESH, MPTHRESH, PTHRESH, INVEPS;
};

class GLvar {
public:
    virtual void devResid(double *ans, const double *mu,
			  const double *pw, const double *y,
			  const int *fac, int n) = 0;
    virtual void varFunc(double *var, const double *mu, int n) = 0;
    virtual const char *distnm() = 0;
};

class GLfamily {
public:
    void initGL(SEXP rho);	/**< initialize from an environment */
    SEXP family;
    const char *fname, *lname;
    GLlink *lnk;
    GLvar *var;
};

class logitlink : public GLlink {
public:
    void linkinv(double *mu, double *muEta, const double* eta, int n);
    const char* nm() {return "logit";}
};

class probitlink : public GLlink {
public:
    void linkinv(double *mu, double *muEta, const double* eta, int n);
    const char* nm() {return "probit";}
};

class clogloglink : public GLlink {
public:
    void linkinv(double *mu, double *muEta, const double* eta, int n);
    const char* nm() {return "cloglog";}
};

class identitylink : public GLlink {
public:
    void linkinv(double *mu, double *muEta, const double* eta, int n);
    const char* nm() {return "identity";}
};

class loglink : public GLlink {
public:
    void linkinv(double *mu, double *muEta, const double* eta, int n);
    const char* nm() {return "log";}
};

class sqrtlink : public GLlink {
public:
    void linkinv(double *mu, double *muEta, const double* eta, int n);
    const char* nm() {return "sqrt";}
};

class constvar : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
    const char *distnm() {return "gaussian";}
};

class mu1muvar : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
    const char *distnm() {return "binomial";}
};

class muvar : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
    const char *distnm() {return "poisson";}
};

class mu2var : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
    const char *distnm() {return "Gamma";}
};

class mu3var : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
    const char *distnm() {return "inverse.gaussian";}
};

#endif /* LME4_GLFAMILY_HPP */

