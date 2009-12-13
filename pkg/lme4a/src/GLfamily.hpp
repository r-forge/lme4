#ifndef LME4_GLFAMILY_HPP
#define LME4_GLFAMILY_HPP

#include "lme4utils.h"		// for definition of SEXP

class GLlink {
public:
    virtual void link(double *etc, const double* mu, int n);
    virtual void linkinv(double *mu, double *muEta,
			 const double* eta, int n);
    static const double LTHRESH, MLTHRESH, MPTHRESH, PTHRESH, INVEPS;
};

class GLvar {
public:
    virtual void devResid(double *ans, const double *mu,
			  const double *pw, const double *y,
			  const int *fac, int n);
    virtual void varFunc(double *var, const double *mu, int n);
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
    void link(double *eta, const double* mu, int n);
    void linkinv(double *mu, double *muEta, const double* eta, int n);
};

class probitlink : public GLlink {
public:
    void link(double *eta, const double* mu, int n);
    void linkinv(double *mu, double *muEta, const double* eta, int n);
};

class clogloglink : public GLlink {
public:
    void link(double *eta, const double* mu, int n);
    void linkinv(double *mu, double *muEta, const double* eta, int n);
};

class identitylink : public GLlink {
public:
    void link(double *eta, const double* mu, int n);
    void linkinv(double *mu, double *muEta, const double* eta, int n);
};

class loglink : public GLlink {
public:
    void link(double *eta, const double* mu, int n);
    void linkinv(double *mu, double *muEta, const double* eta, int n);
};

class sqrtlink : public GLlink {
public:
    void link(double *eta, const double* mu, int n);
    void linkinv(double *mu, double *muEta, const double* eta, int n);
};

class constvar : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
};

class mu1muvar : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
};

class muvar : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
};

class mu2var : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
};

class mu3var : public GLvar {
public:
    void devResid(double *ans, const double *mu,
		  const double *pw, const double *y,
		  const int *fac, int n);
    void varFunc(double *var, const double *mu, int n);
};

#endif /* LME4_GLFAMILY_HPP */

