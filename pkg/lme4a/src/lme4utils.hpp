#ifndef LME4_LME4UTILS_HPP
#define LME4_LME4UTILS_HPP

/** Non-inlined utilities
 */

SEXP CHM_SP2SEXP(CHM_SP A, const char *cls);
SEXP CHM_SP2SEXP(CHM_SP A, const char *cls, const char *uplo);
SEXP CHM_SP2SEXP(CHM_SP A, const char *cls, const char *uplo, const char *diag);
double *VAR_REAL_NULL(SEXP rho, SEXP nm, int len, int nullOK,
		      int absentOK);
double *VAR_REAL_NULL(SEXP rho, SEXP nm, int len, int nullOK);
double *VAR_REAL_NULL(SEXP rho, SEXP nm, int len);
CHM_SP VAR_CHM_SP(SEXP rho, SEXP nm, int nrow, int ncol);
CHM_FR VAR_CHM_FR(SEXP rho, SEXP nm, int n);
double *VAR_dMatrix_x(SEXP rho, SEXP nm, int nrow, int ncol);

CHM_SP CHM_SP_copy_in_place(CHM_SP dest, CHM_SP src);

SEXP getListElement(SEXP list, SEXP names, const char *str);

#endif /* LME4_LME4UTILS_HPP */
