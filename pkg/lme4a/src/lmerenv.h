#ifndef LME4_LMERENV_H
#define LME4_LMERENV_H

#include <R.h>
// Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
// GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc.
#include <Rdefines.h>

#ifdef	__cplusplus
#define NO_C_HEADERS
#include <cstdio>
#include <cmath>
#include <cstring>
#include "merenv.h"
#include "Matrix.h"

class lmerenv : public merenv {
public:
    lmerenv(SEXP rho);		//< instantiate from an environment
//    ~lmerenv(){}
    double update_dev(SEXP thnew);
    int validate();

private:
    int REML;
    double *RX, *RZX, *XtX, *Xty, *ZtX, *Zty, *ldRX2;
};
#endif /* __cplusplus */

#endif /* LME4_LMERENV_H */
