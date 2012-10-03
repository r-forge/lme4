############           
## FIXME: use xapply?
results <- vector("list",
                  length(rhobegvec)*length(begendvec)*length(tolPwrssvec))
  ctr <- 0
  for (i in seq_along(rhobegvec)) {
    for (j in seq_along(begendvec)) {
      for (k in seq_along(tolPwrssvec)) {
        ctr <- ctr+1
        cat("*", i,j,k,ctr,"\n")
        cur.opt <- glmer(form,Contraception,binomial,
                         tolPwrss=tolPwrssvec[k],
                         control=list(rhobeg=rhobegvec[i],
                           rhoend=rhobegvec[i]*begendvec[j]))
        if (inherits(cur.opt,"try-error")) {
          d <- dbad
          d$result <- attr(cur.opt,"condition")$message
        } else {
          d <- d0
          d <- c(d,list(options=c(control,list(tolPwrss=tolPwrssvec[k])),
                        time=sum(t0[1:2]), ## user +system
                        parameters=allcoef(cur.opt),
                        deviance=deviance(cur.opt),
                        KKT=c(NA,NA),
                        ## FIXME: extract KKT from lme4Eigen?
                        bad=NA,
                        result=NA_character_))
        }
        results[[ctr]] <- d
        save("results","sessinfo",file=fn)
      }
    }
  }
