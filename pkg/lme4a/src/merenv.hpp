#ifndef LME4_MERENV_H
#define LME4_MERENV_H

#ifdef	__cplusplus
#define NO_C_HEADERS
#include <cstdio>
#include <cmath>
#include <cstring>
extern "C" {
#endif 

#include <R.h>
// Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
// GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc.
#include <Rdefines.h>
#include "Matrix.h"

    SEXP glmerenv_linkinv(SEXP rho);

    SEXP lme4_dup_env_contents(SEXP dest, SEXP src, SEXP nms);

    SEXP lmerenv_deviance(SEXP rho, SEXP newth);
    SEXP lmerenv_validate(SEXP rho);

    SEXP merenvtrms_condVar(SEXP rho, SEXP scale);
    SEXP merenvtrms_show(SEXP rho);
    SEXP merenvtrms_validate(SEXP rho);

#ifdef	__cplusplus
}

static int i1 = 1;
static double one = 1, mone = -1, zero = 0;

class merenv {
    /// Basic class for mixed-effects.
public:
    ~merenv(){
	delete L;
	delete Lambda;
	delete Ut;
	delete Zt;
    }
/** 
 * Initialize from an environment.
 * 
 * @param rho pointer to an environment
 */
    void initMer(SEXP rho);	
/** 
 * Update the linear predictor.
 * 
 * Update the linear predictor using the offset, if present, and the
 * random-effects.  Updates from the fixed-effects are done in the
 * derived classes. 
 */
    void update_eta_Ut();

/** 
 *  Update Lambda and Ut from theta.
 * 
 * @param thnew pointer to a numeric vector or new values for theta
 */
    void update_Lambda_Ut(SEXP thnew);
/** 
 * Update the penalized residual sum-of-squares.
 *
 * @return updated penalized residual sum-of-squares.
 */
    double update_prss();
/** 
 * Create the crossproduct of Lambda and src in ans.
 * 
 * @param src dense matrix with q rows
 * @param ans matrix of same size as src to be overwritten
 * 
 * @return ans, overwritten with the product
 */
    CHM_DN crossprod_Lambda(CHM_DN src, CHM_DN ans);
/** 
 * Return the crossproduct of Lambda and src
 * 
 * @param src sparse matrix with q rows
 * 
 * @return crossprod(Lambda, src)
 */
    CHM_SP spcrossprod_Lambda(CHM_SP src);
/** 
 * Solve L ans = P src (dense case)
 * 
 * @param src dense matrix of values on the right hand side
 * 
 * @return a dense matrix of the same size as src
 */
    CHM_DN solvePL(CHM_DN src);
/** 
 * Solve L ans = P src (sparse case)
 * 
 * @param src sparse matrix of values on the right hand side
 * 
 * @return a sparse matrix of the same size as src
 */
    CHM_SP solvePL(CHM_SP src);

    int
	N,		/**< length of eta (can be a multiple of n) */
	n,		/**< number of observations */
	p,		/**< number of fixed effects */
	q;		/**< number of random effects */
    double
	*Lambdax,	     /**< x slot of Lambda */
	*eta,		     /**< linear predictor values */
	*fixef,		     /**< fixed-effects parameters */
	*ldL2,		     /**< log-determinant of L squared */
	*prss,		     /**< penalized residual sum-of-squares */
	*theta,		     /**< parameters that determine Lambda */
	*u,		     /**< unit random-effects vector */
	*weights,	     /**< prior weights (may be NULL) */
	*y;		     /**< response vector */
    CHM_FR
	L;			/**< sparse Cholesky factor */
    CHM_SP
	Lambda,	     /**< relative covariance factor */
	Ut,	     /**< model matrix for unit random effects */
	Zt;	     /**< model matrix for random effects */
    int *Lind,	     /**< Lambdax index vector into theta (1-based) */
	nLind;	     /**< length of Lind */

private:
    int nth;	      /**< number of elements of theta */
    double *offset;   /**< offset of linear predictor (may be NULL) */
};

class mersparse : virtual public merenv { // merenv with sparse X
public:
    ~mersparse() {
	delete X;
	delete RX;
	delete RZX;
    }
    void initMersd(SEXP rho);	/**< initialize from environment */
    void update_eta();		/**< update linear predictor */
    CHM_SP X,			/**< model matrix for fixed effects */
	RX,		     /**< Cholesky factor for fixed-effects */
	RZX;		     /**< cross-product in Cholesky factor */
};

class merdense : virtual public merenv { // merenv with dense X
public:
    void initMersd(SEXP rho);	/**< initialize from environment */
    void update_eta();		/**< update linear predictor */
    double 
	*RX,	       /**< upper Cholesky factor for fixed-effects */
	*RZX,	       /**< cross-product in Cholesky factor */
	*X;	       /**< model matrix for fixed effects */
};

class lmer : virtual public merenv { // components common to LMMs
public:
    void initLMM(SEXP rho);	/**< initialize from environment */
    void LMMdev1();   /**< initial shared part of deviance update */
    void LMMdev2();   /**< secondary shared part of deviance update */
    double LMMdev3(); /**< tertiary shared part of deviance update */
    int
	REML;			/**< logical - use REML? */
    double
	*Xty,			/**< cross product of X and y */
	*Zty,			/**< cross product of Z and y */
	*ldRX2;			/**< log-determinant of RX squared */
    CHM_DN
	cu;		  /**< Intermediate value in solution for u */
};

class lmerdense : public merdense, public lmer {
public:
    lmerdense(SEXP rho);	//< construct from an environment
    double update_dev(SEXP thnew);
    int validate() {
	return 1; 
    }
    double *XtX, *ZtX;
};

class lmersparse : public mersparse, public lmer {
public:
    lmersparse(SEXP rho);	//< construct from an environment
    ~lmersparse(){
	delete XtX;
	delete ZtX;
    }
    double update_dev(SEXP thnew);
    int validate() {
	return 1;
    }
    CHM_SP XtX, ZtX;
};

class merenvtrms : public merenv {
public:
    merenvtrms(SEXP rho);  	//< construct from an environment
    
    SEXP condVar(double scale);	//< create the conditional variance array
    void show();		//< debugging output

    int validate() {
	return 1;
    }
    SEXP flist;			//< pointer to list of grouping factors
    int nfac;			//< number of grouping factors
    int *nl;			//< number of levels per factor
    int *apt;			//< rle pointers into assign array
    int ntrm;			//< number of terms
    int *nc;			//< number of columns per term
};

#include "GLfamily.hpp"

static logitlink logitlnk;
static probitlink probitlnk;
static identitylink identitylnk;
static loglink loglnk;
static sqrtlink sqrtlnk;
static constvar constvr;
static mu1muvar mu1muvr;
static muvar muvr;
static mu2var mu2vr;
static mu3var mu3vr;

class glmerenv : virtual public merenv { // components common to GLMMs
public:
    glmerenv(SEXP rho);		//< initialize from an environment
    SEXP family;
    double
	*mu,			//< conditional means of response
	*muEta,			//< diagonal of d mu/d eta
	*sqrtrwt,		//< square root of resid weights
	*var,			//< conditional variances of response
	devres;
    GLlink *ll;
    GLvar *vv;
    void linkinv() {ll->linkinv(mu, muEta, eta, n);}
    void devResid() {vv->devResid(&devres, mu, weights, y,
				  (int*)0, n);}
    void varFunc() {vv->varFunc(var, mu, n);}
    const char *fam;
    const char *lnk;
};


#endif /* __cplusplus */

#endif /* LME4_MERENV_H */

