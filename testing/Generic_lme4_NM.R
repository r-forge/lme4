source("xapply.R")
source("miscfuns.R")
require(plyr)

lme4fit <- function(Problem,Data,Form,Family,
                     method="NM",
                     tolPwrssvec=10^seq(-11,-5,by=0.25),
                     rhobegvec=2*10^(-8:0),
                     begendvec=10^(-3:-6),
                     nAGQvec=c(2^(0:4),25),...)  ## lme4 has max nAGQ=25
{
  require(lme4)
  if (method=="NM") {
    results <- xapply(tolPwrssvec,nAGQvec,FUN=dfun_lme4,
                      MoreArgs=list(problem=Problem,
                        form=Form,data=Data,family=Family),
                      .progress="text",pbargs=list(style=3))
  } else {
    results <- xapply(FUN=dfun_lme4Eigen_bobyqa,
                      rhobegvec,begendvec,tolPwrssvec,nAGQvec,
                      MoreArgs=list(problem=Problem,
                        form=Form,data=Data,family=Family),
                      .progress="txt",pbargs=list(style=3))
  }
  detach("package:lme4")
  results
}
