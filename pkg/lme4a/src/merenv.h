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

SEXP lmerenv_deviance(SEXP rho, SEXP newth);
SEXP lmerenv_validate(SEXP rho);

#ifdef	__cplusplus
}

static int i1 = 1;
static double one = 1, mone = -1, zero = 0;

class merenv {
public:
    ~merenv(){
	delete L;
	delete Lambda;
	delete Ut;
	delete Zt;
    }
    void initMer(SEXP rho);
    void update_eta_Ut();
    void update_Lambda_Ut(SEXP thnew);
    double update_prss();
    CHM_DN crossprod_Lambda(CHM_DN src, CHM_DN ans);
    CHM_SP spcrossprod_Lambda(CHM_SP src);
    CHM_DN solvePL(CHM_DN src);
    CHM_SP solvePL(CHM_SP src);
    int N, n, p, q;
    double *Lambdax, *eta, *fixef, *ldL2, *prss, *theta, *u, *weights, *y;
    CHM_FR L;
    CHM_SP Lambda, Ut, Zt;
    int *Lind, nLind;

private:
    int nth;
    double *offset;
};

class mersparse : virtual public merenv { // merenv with sparse X
public:
    ~mersparse() {
	delete X;
	delete RX;
	delete RZX;
    }
    void initMersd(SEXP rho);
    void update_eta();
    CHM_SP X, RX, RZX;
};

class merdense : virtual public merenv { // merenv with dense X
public:
    void initMersd(SEXP rho);
    void update_eta();
    double *X, *RX, *RZX;
};

class lmer : virtual public merenv { // components common to LMMs
public:
    void initLMM(SEXP rho);
    void LMMdev1();
    void LMMdev2();
    double LMMdev3();
    int REML;
    double *Xty, *Zty, *ldRX2;
    CHM_DN cu;
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

#endif /* __cplusplus */

#endif /* LME4_MERENV_H */
