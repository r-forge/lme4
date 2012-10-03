source("xapply.R")
source("miscfuns.R")
require(plyr)

## lme4a does not do AGQ ...
lme4afit <- function(Problem,Data,Form,Family) 
{
  require(lme4a)
  results <- list(dfun_lme4other(nAGQ=1,
                                 problem=Problem,
                                 pkg="lme4a",
                                 form=Form,data=Data,family=Family))
  detach("package:lme4a")
  results
}
