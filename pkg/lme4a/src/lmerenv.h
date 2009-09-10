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

class lmerenv : public merenv {
public:
    lmerenv(SEXP rho);		//< instantiate from an environment
    ~lmerenv(){
	delete sX;
    }
    void update_dev(SEXP thnew);
    int validate();

private:
    double *RX, *RZX, *XtX, *Xty, *ZtX, *Zty;
};
#endif /* __cplusplus */

#endif /* LME4_LMERENV_H */
