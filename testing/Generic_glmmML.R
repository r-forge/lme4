source("xapply.R")
source("miscfuns.R")
require(plyr)

glmmMLfit <- function(Problem,Data,Form,Family,nAGQvec=2^(0:5)) {
  require(glmmML)
  ## HACK!! glmmML evaluates the call in the environment of the
  ##   formula, so we have to stuff the information in there ...
  assign("Data",Data,environment(Form))
  assign("Form",Form,environment(Form))
  assign("cluster",Data$cluster,environment(Form))
  g1 <-  glmmML(Form, family=Family, data=Data, cluster = cluster)
    
  t1 <- system.time(g1 <- glmmML(Form, family=Family, data=Data,
                                 cluster = cluster))["elapsed"]
  r1  <- list(list(problem=Problem,
                   pkg="glmmML",
                   method="laplace",
                   nAGQ=1,options=NULL,
                   time=t1,  ## use ELAPSED time
                   parameters=allcoef(g1),
                   deviance=deviance(g1),
                   options=NULL))

  dfun <- function(a) {
    t1 <- system.time(g1 <- glmmML(Form, family=Family, data=Data,
                                 cluster = cluster,method="ghq",
                                   n.points=a))["elapsed"]
    list(problem=Problem,
         pkg="glmmML",
       method="ghq",
       nAGQ=a,options=NULL,
         time=t1,  ## use ELAPSED time
         parameters=allcoef(g1),
         deviance=deviance(g1),
         options=NULL)
  }

  r3 <- llply(as.list(nAGQvec),
              .progress="text",
              .fun = dfun)
  
  results <- c(r1,r3)
  return(results)
}
save("results",file=fn)
