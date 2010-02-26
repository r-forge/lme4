#ifndef LME4_MERENV_HPP
#define LME4_MERENV_HPP
#include "lme4utils.h"		// to get definitions of SEXP
#include "matrix.hpp"

class merenv {
public:
    merenv(SEXP rho);
    ~merenv();
/** 
 * Update the linear predictor.
 * 
 * Update the linear predictor using the offset, if present, the fixed
 * effects and the random effects.
 */
    void update_gamma();
/** 
 *  Update Lambda, and Ut from theta and sqrtXwt
 * 
 * @param thnew pointer to a numeric vector or new values for theta
 */
    void update_Lambda_Ut(SEXP thnew);
/** 
 * Update the penalized residual sum-of-squares.
 *
 * @return updated penalized residual sum-of-squares.
 */
    double update_pwrss();
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
    CHM_r *Vp(); /**< return pointer to a copy of X with rows scaled by sqrtXwt */

    int
	N,		/**< length of gamma (can be a multiple of n) */
	n,		/**< number of observations */
	nth,		/**< number of elements of theta */
	p,		/**< number of fixed effects */
	q,		/**< number of random effects */
	diagonalLambda,
	sparseX,
	verbose;
    
    double
	*fixef,	      /**< fixed-effects parameters */
	*gam,	      /**< linear predictor (not called gamma b/c
			 of conflicts with the gamma function) */ 
	*Lambdax,     /**< x slot of Lambda matrix */
	*ldL2,	      /**< log-determinant of L squared */
	*mu,	      /**< conditional mean response */
	*offset,      /**< offset of linear predictor (may be NULL) */
	*pwrss,	      /**< penalized, weighted RSS */
	*sqrtrwt,     /**< sqrt of residual weights */
	*sqrtXwt,     /**< sqrt of weights for X->V and U->L */  
	*theta,	      /**< parameters that determine Lambda */
	*u,	      /**< unit random-effects vector */
	*weights,     /**< prior weights (may be NULL) */
	*y;	      /**< response vector */
    CHM_FR L;	      /**< sparse Cholesky factor */
    CHM_SP Lambda,    /**< relative covariance factor */
	Ut,	      /**< model matrix for unit random effects */
	Zt;	      /**< model matrix for random effects */
    int *Lind,	      /**< lambdax index vector into theta (1-based) */
	nLind;	      /**< length of Lind */
    CHM_r *Xp, *RZXp;
    Cholesky_r *RXp;
};

/// Linear mixed-effects model
class lmer : public merenv {
public:
    lmer(SEXP);   /**< initialize from environment */
    ~lmer();
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

/// generalized linear mixed models

class glmer : public merenv {
public:
    glmer(SEXP rho);		/**< initialize from an environment */
    SEXP family;
    double 
	*muEta,			/**< diagonal of d mu/d eta */
	*var,			/**< conditional variances of response */
	*wtres,			/**< weighted residuals */
	devres;
    GLfamily fam;               /**< GLM family of functions */
    double PIRLS();		/**< deviance at updated u */
    double PIRLSbeta();		/**< deviance at updated u and beta */
    double IRLS();		/**< deviance at updated beta */
    double Laplace();		/**< Laplace approximation to the deviance */
    void update_sqrtrwt();
    void update_sqrtXwt();
    double update_wtres();
    void link() {fam.lnk->link(gam, mu, n);}
    void linkinv() {fam.lnk->linkinv(mu, muEta, gam, n);}
    void devResid() {fam.var->devResid(&devres, mu, weights, y,
				      (int*)0, n);}
    void varFunc() {fam.var->varFunc(var, mu, n);}
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

#endif /* LME4_MERENV_HPP */
