##
## substitute DATNAME before NAME, FIXFORM/RANFORM before FORM, ..
s/DATNAME/cbpp/g
s/FIXNAMES/fixef..Intercept.,fixef.period2,fixef.period3,fixef.period4/g
s/NAME/Cbpp/g
s/DATPKG/lme4/g
s/FAMILY/binomial/g
s/FIXFORM/cbind(incidence,size-incidence)~period/g
s/RANFORM/1|herd/g
s/FORM/cbind(incidence,size-incidence)~period+(1|herd)/g
s/CLUSTER/herd/g
s/RESPONSE/incidence/g

