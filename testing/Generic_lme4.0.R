source("xapply.R")
source("miscfuns.R")
require(plyr)

lme4.0fit <- function(Problem,Data,Form,Family,
                     nAGQvec=2^(0:5))
{
  require(lme4.0)
  results <- xapply(nAGQvec,FUN=dfun_lme4other,
                    MoreArgs=list(problem=Problem,
                      pkg="lme4.0",
                      form=Form,data=Data,family=Family),
                      .progress="text",pbargs=list(style=3))
  detach("package:lme4.0")
  results
}
