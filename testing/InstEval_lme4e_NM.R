fn <- "InstEval_lme4e_NM.RData"
source("Generic_lme4e.R")

sessinfo <- sessionInfo()

data(InstEval,package="lme4Eigen")
results <- lme4efit("InstEval",
                    Data=InstEval,
                    Form=y ~ 1 + (1|s)+(1|d) + (1|dept:service),
                    method="NM",
                    nAGQvec=NA, ## turn off GLMM stuff
                    tolPwrssvec=NA,
                    REML=FALSE)
save("results",file=fn)
