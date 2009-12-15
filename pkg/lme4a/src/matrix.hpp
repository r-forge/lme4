#include "lme4utils.h"		// for CHM_DN, CHM_SP, c

// Classes of double precision matrices

/** abstract class of double precision matrices in CHOLMOD form */
class CHM_r {
public:
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y) = 0;
};

/** concrete class of double precision dense matrices */
class CHM_rd : public CHM_r {
public:
    CHM_rd(SEXP x);
    ~CHM_rd(){delete A;}
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y);
protected:
    CHM_DN A;
};

/** concrete class of double precision general matrices */
class CHM_rs : public CHM_r {
public:
    CHM_rs(SEXP x);
    ~CHM_rs(){delete A;}
    virtual void drmult(int transpose, double alpha, double beta,
			CHM_DN X, CHM_DN Y);
protected:
    CHM_SP A;
};
