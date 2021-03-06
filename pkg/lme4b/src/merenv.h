#ifndef LME4_MERENV_H
#define LME4_MERENV_H
#include "lme4utils.h"

#ifdef __cplusplus
extern "C" {
#endif

    SEXP glmer_IRLS(SEXP rho);
    SEXP glmer_PIRLS(SEXP rho);
    SEXP glmer_PIRLSbeta(SEXP rho);
    SEXP glmer_update_RX(SEXP rho);
    
    SEXP lme4_dup_env_contents(SEXP dest, SEXP src, SEXP nms);

    SEXP lmer_deviance(SEXP rho, SEXP newth);
    SEXP lmer_validate(SEXP rho);

    SEXP lmerold_validate(SEXP rho);
    SEXP lmerold_deviance(SEXP rho, SEXP newth);

    SEXP merenv_update_Lambda_Ut(SEXP rho, SEXP newth);

    SEXP merenvtrms_condVar(SEXP rho, SEXP scale);
    SEXP merenvtrms_show(SEXP rho);
    SEXP merenvtrms_validate(SEXP rho);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* LME4_MERENV_H */
