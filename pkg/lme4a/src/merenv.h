#ifndef LME4_MERENV_H
#define LME4_MERENV_H

#ifdef	__cplusplus
#define NO_C_HEADERS
#include <cstdio>
#include <cmath>
#include <cstring>
#endif 

#include <R.h>
// Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
// GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc.
#include <Rdefines.h>
#include "Matrix.h"

#ifdef	__cplusplus

class merenv {
public:
    merenv(SEXP rho);		//< construct from an environment
    ~merenv(){
	delete L;
	delete Lambda;
	delete Ut;
	delete Zt;
    }
    static int i1;
    int N, n, p, q;
    double *Lambdax, *eta, *fixef, *ldL2, *prss, *theta, *u, *weights, *y;
    CHM_FR L;
    CHM_SP Lambda, Ut, Zt;
    int *Lind, nLind;

private:
    int nth;
    double *offset;
};

class mersparse : public merenv { // merenv with sparse X
public:
    mersparse(SEXP rho);	//< construct from an environment
    ~mersparse() {
	delete X;
	delete RX;
	delete RZX;
    }
    void update_eta();
    CHM_SP X, RX, RZX;
};

class merdense : public merenv { // merenv with dense X
public:
    merdense(SEXP rho);		//< construct from an environment
    void update_eta();
    double *X, *RX, *RZX;
};

class lmer {			// components common to LMMs
public:
    lmer(SEXP rho);
    int REML;
    double *Xty, *Zty, *ldRX2;
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
