#ifndef LME4_MERENV_HPP
#define LME4_MERENV_HPP
#include "lme4utils.h"		// to get definitions of SEXP
#include "matrix.hpp"

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
	delete Xp;
	delete RXp;
	delete RZXp;
    }
    void initMer(SEXP rho);

/** 
 * Update the linear predictor.
 * 
 * Update the linear predictor using the offset, if present, and the
 * random-effects.  Updates from the fixed-effects are done in the
 * derived classes. 
 */
    void update_gamma();

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
 * @return ans, overwritten with the product
 */
    CHM_DN crossprod_Lambda(CHM_DN src, CHM_DN ans);
/** 
 * Return the crossproduct of Lambda and src
 * 
 * @param src sparse matrix with q rows
 * @return crossprod(Lambda, src)
 */
    CHM_SP spcrossprod_Lambda(CHM_SP src);
/** 
 * Solve L ans = P src (dense case)
 * 
 * @param src dense matrix of values on the right hand side
 * @return a dense matrix of the same size as src
 */
    CHM_DN solvePL(CHM_DN src);
/** 
 * Solve L ans = P src (sparse case)
 * 
 * @param src sparse matrix of values on the right hand side
 * @return a sparse matrix of the same size as src
 */
    CHM_SP solvePL(CHM_SP src);

    int
	N,		/**< length of gamma (can be a multiple of n) */
	n,		/**< number of observations */
	p,		/**< number of fixed effects */
	q,		/**< number of random effects */
	sparseX;
    double
	*Lambdax,	     /**< x slot of Lambda */
	*fixef,		     /**< fixed-effects parameters */
	*gam,		     /**< linear predictor values (not called
				gamma b/c that conflicts with the gamma func */
	*ldL2,		     /**< log-determinant of L squared */
	*mu,		     /**< conditional mean response */
	*prss,		     /**< penalized residual sum-of-squares */
	*sqrtrwt,	     /**< square root of residual weights */
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
    CHM_r *Xp, *RZXp;
    Cholesky_r *RXp;

private:
    int nth;	      /**< number of elements of theta */
    double *offset;   /**< offset of linear predictor (may be NULL) */
};

class mersparse : virtual public merenv { // merenv with sparse X
public:
    mersparse(SEXP rho);
    ~mersparse() {
	delete X;
	delete RX;
	delete RZX;
    }
    CHM_SP X,		     /**< model matrix for fixed effects */
	RZX;		     /**< cross-product in Cholesky factor */
    CHM_FR RX;		     /**< Cholesky factor for fixed-effects */
};

// merdense and mersparse are no longer needed.  Pull when lmerdense
// and lmersparse are removed.
class merdense : virtual public merenv { // merenv with dense X
public:
    merdense(SEXP rho);
    double 
	*RX,	       /**< upper Cholesky factor for fixed-effects */
	*RZX,	       /**< cross-product in Cholesky factor */
	*X;	       /**< model matrix for fixed effects */
};

class lmer : virtual public merenv { // components common to LMMs
public:
    lmer(SEXP rho);   /**< initialize from environment */
    ~lmer() {
	delete XtXp;
	delete ZtXp;
    }
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
    dMatrix *XtXp;
    CHM_r *ZtXp;
};

class lmerdense : public lmer, public merdense {
public:
    lmerdense(SEXP rho);	/**< construct from an environment */
    double update_dev(SEXP thnew);
    int validate() {return 1;}
    double *XtX, *ZtX;
};

class lmersparse : public lmer, public mersparse {
public:
    lmersparse(SEXP rho);	/**< construct from an environment */
    ~lmersparse(){
	delete XtX;
	delete ZtX;
    }
    double update_dev(SEXP thnew);
    int validate() {return 1;}
    CHM_SP XtX, ZtX;
};

class merenvtrms : public merenv {
public:
    merenvtrms(SEXP rho);  	/**< construct from an environment */
    
    SEXP condVar(double scale);	/**< create the conditional variance array */
    void show();		/**< debugging output */

    int validate() {
	return 1;
    }
    SEXP flist;			/**< pointer to list of grouping factors */
    int nfac;			/**< number of grouping factors */
    int *nl;			/**< number of levels per factor */
    int *apt;			/**< rle pointers into assign array */
    int ntrm;			/**< number of terms */
    int *nc;			/**< number of columns per term */
};

class mernew {
public:
    mernew(SEXP rho);
    ~mernew();
/** 
 * Update the linear predictor.
 * 
 * Update the linear predictor using the offset, if present, the fixed
 * effects and the random effects.
 */
    void update_gamma();
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
 * @return ans, overwritten with the product
 */
    CHM_DN crossprod_Lambda(CHM_DN src, CHM_DN ans);
/** 
 * Return the crossproduct of Lambda and src
 * 
 * @param src sparse matrix with q rows
 * @return crossprod(Lambda, src)
 */
    CHM_SP spcrossprod_Lambda(CHM_SP src);
/** 
 * Solve L ans = P src (dense case)
 * 
 * @param src dense matrix of values on the right hand side
 * @return a dense matrix of the same size as src
 */
    CHM_DN solvePL(CHM_DN src);
/** 
 * Solve L ans = P src (sparse case)
 * 
 * @param src sparse matrix of values on the right hand side
 * @return a sparse matrix of the same size as src
 */
    CHM_SP solvePL(CHM_SP src);

    int
	N,		/**< length of gamma (can be a multiple of n) */
	n,		/**< number of observations */
	p,		/**< number of fixed effects */
	q,		/**< number of random effects */
	sparseX;
    double
	*Lambdax,	     /**< x slot of Lambda */
	*fixef,		     /**< fixed-effects parameters */
	*gam,		     /**< linear predictor (not called gamma b/c
				that conflicts with the gamma function) */
	*ldL2,		     /**< log-determinant of L squared */
	*mu,		     /**< conditional mean response */
	*prss,		     /**< penalized residual sum-of-squares */
	*sqrtrwt,	     /**< sqrt of residual weights */
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
    CHM_r *Xp, *RZXp;
    Cholesky_r *RXp;

private:
    int nth;	      /**< number of elements of theta */
    double *offset;   /**< offset of linear predictor (may be NULL) */
};

/// Linear mixed-effects model
class lmernew : public mernew {
public:
    lmernew(SEXP);   /**< initialize from environment */
    ~lmernew() {
	delete XtXp;
	delete ZtXp;
    }
    double update_dev(SEXP thnew);
    int validate(){return 1;}
    int
	REML;			/**< logical - use REML? */
    double
	*Xty,			/**< cross product of X and y */
	*Zty,			/**< cross product of Z and y */
	*ldRX2;			/**< log-determinant of RX squared */
    CHM_DN
	cu;		  /**< Intermediate value in solution for u */
    dMatrix *XtXp;
    CHM_r *ZtXp;
};

#include "GLfamily.hpp"
class glmer : virtual public merenv { // components common to GLMMs
public:
    glmer(SEXP rho);		/**< initialize from an environment */
    SEXP family;
    double 
	*eta,			/**< conditional mean on link scale */
	*muEta,			/**< diagonal of d mu/d eta */
	*var,			/**< conditional variances of response */
	*wtres,			/**< weighted residuals at current estimates */
	devres;
    GLfamily fam;
    double PIRLS();		/**< deviance at updated u */
    double PIRLSbeta();		/**< deviance at updated u and beta */
    double IRLS();		/**< deviance at updated beta */
    void update_sqrtrwt();
    void link() {fam.lnk->link(eta, mu, n);}
    void linkinv() {fam.lnk->linkinv(mu, muEta, eta, n);}
    void devResid() {fam.var->devResid(&devres, mu, weights, y,
				      (int*)0, n);}
    void varFunc() {fam.var->varFunc(var, mu, n);}
};

class glmerdense : public glmer, public merdense {
public:
    glmerdense(SEXP rho);	/**< construct from an environment */
    ~glmerdense() {delete[] V;}

    void update_V();
    double update_dev(SEXP thnew);
    int validate() {
	return 1; 
    }
    double *V;
    double *XtX, *ZtX;
};

class glmersparse : public glmer, public mersparse {
public:
    glmersparse(SEXP rho);	/**< construct from an environment */
    ~glmersparse(){
	delete XtX;
	delete ZtX;
	M_cholmod_free_sparse(&V, &c);
	delete V;
    }
    void update_V();
    double update_dev(SEXP thnew);
    int validate() {
	return 1;
    }
    CHM_SP XtX, ZtX, V;
};
    
#endif /* LME4_MERENV_HPP */
