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
    merenv(SEXP rho);		//< instantiate from an environment
    ~merenv(){
	delete L;
	delete Lambda;
	delete Ut;
	delete Zt;
	delete sX;
    }
    void update_Lambda_Ut(SEXP thnew);
    static int i1;
    int validate() {		// validation occurs in constructor
	return 1;
    }
    CHM_DN crossprod_Lambda(CHM_DN rhs, CHM_DN ans);
    void update_eta();
    int N, n, p, q;
    double *X, *eta, *fixef, *ldL2, *prss, *theta, *u, *weights, *y;
    CHM_FR L;
    CHM_SP Lambda, Ut, Zt, sX;

private:
    int *Lind, nLind, nth;
    double *Lambdax, *offset;
};
#endif /* __cplusplus */

#endif /* LME4_MERENV_H */
