##
## substitute DATNAME before NAME, FIXFORM/RANFORM before FORM, ..
s/DATNAME/bacteria/g
s/FIXNAMES/fixef..Intercept.,fixef.trtdrug,fixef.trtdrug.,fixef.I.week...2.TRUE/g
s/NAME/Bacteria/g
s/DATPKG/MASS/g
s/FIXFORM/y~trt+I(week>2)/g
s/RANFORM/1|ID/g
s/FORM/y~trt+I(week>2)+(1|ID)/g
s/FAMILY/"binomial"/g
s/CLUSTER/ID/g
s/RESPONSE/y/g

