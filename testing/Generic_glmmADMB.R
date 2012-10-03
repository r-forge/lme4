## TO DO: add importance sampling, GHQ where possible?

source("xapply.R")
source("miscfuns.R")
require(plyr)

glmmADMBfit <- function(Problem,Data,Form,Family,...) 
{
  require(glmmADMB)
  if (is.function(Family)) Family <- deparse(substitute(Family))
  gt0 <- system.time(g0 <- glmmadmb(formula=Form,
                                 data=Data,family=Family,...))

  results <- list(list(problem=Problem,
     pkg="glmmADMB",
     method="ADMB",
     nAGQ=1,options=NULL,
     time=gt0[3],  ## use ELAPSED time
     parameters=allcoef(g0),
     deviance=-2*c(logLik(g0)),
     options=NULL))

  detach("package:glmmADMB")
  results
}
