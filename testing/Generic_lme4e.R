source("xapply.R")
source("miscfuns.R")
require(plyr)

lme4efit <- function(Problem,Data,Form,Family=gaussian,
                     method="NM",
                     tolPwrssvec=10^seq(-11,-5,by=0.25),
                     rhobegvec=2*10^(-8:0),
                     begendvec=10^(-3:-6),
                     nAGQvec=c(2^(0:4),25),  ## lme4Eigen has max nAGQ=25
                     ...)
{
  require(lme4Eigen)
  if (method=="NM") {
    results <- xapply(tolPwrssvec,nAGQvec,FUN=dfun_lme4Eigen_NM,
                      MoreArgs=list(problem=Problem,
                        form=Form,data=Data,family=Family,...),
                      .progress="text",pbargs=list(style=3))
  } else {
    results <- xapply(FUN=dfun_lme4Eigen_bobyqa,
                      rhobegvec,begendvec,tolPwrssvec,nAGQvec,
                      MoreArgs=list(problem=Problem,
                        form=Form,data=Data,family=Family,...),
                      .progress="txt",pbargs=list(style=3))
  }
  detach("package:lme4Eigen")
  results
}
