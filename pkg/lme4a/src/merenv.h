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
extern "C" {
#endif

    SEXP lmerenv_deviance(SEXP rho, SEXP thnew);
    SEXP lmerenv_validate(SEXP rho);

#ifdef	__cplusplus
}

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
    int validate();

private:
    static int i1;
    int *Lind, N, n, nLind, nth, p, q;
    double *Lambdax, *X, *beta, *eta, *ldL2, *offset,
	*prss, *theta, *u, *weights, *y;
    CHM_FR L;
    CHM_SP Lambda, Ut, Zt, sX;
    void update_eta();
    CHM_DN crossprod_Lambda(CHM_DN rhs, CHM_DN ans);
};
#endif /* __cplusplus */

#endif /* LME4_MERENV_H */
