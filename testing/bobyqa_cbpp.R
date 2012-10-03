fn <- "bobyqa_cbpp.RData"

library(lme4Eigen)
library(optimx)
sessinfo <- sessionInfo()

data(cbpp)
Data <- cbpp
Family <- binomial

form <- cbind(incidence, size - incidence) ~ period + (1 | herd)

## extract the objective function for the second stage optimization
ff <- glmer(form, Data, Family, devFunOnly=2L)

val <- c(get("pp", environment(ff))$theta, get("pp", environment(ff))$beta(0))

d0 <- list(problem="cbpp",method="bobyqa",optimx=TRUE)
dbad <- c(d0,as.list(c(options=NA,time=NA,parameters=NA,deviance=NA,KKT=NA,bad=NA,result=NA)))
          
rhobegvec <- 2*10^(-8:0)
begendvec <- 10^(-3:-6)
tolPwrssvec <- 10^seq(-9,-5)

results <- vector("list",
                  length(rhobegvec)*length(begendvec)*length(tolPwrssvec))
  ctr <- 0
  for (i in seq_along(rhobegvec)) {
    for (j in seq_along(begendvec)) {
      for (k in seq_along(tolPwrssvec)) {
        ## reset ff each time ... just in case it gets corrupted
        ff <- glmer(form, Data, Family, devFunOnly=2L)
        ctr <- ctr+1
        cat("*", i,j,k,ctr,"\n")
        control <- list(rhobeg=rhobegvec[i],rhoend=rhobegvec[i]*begendvec[j])
        ## tolPwrss is read from ff's environment
        tolPwrss0 <- get("tolPwrss",envir=environment(ff))
        assign("tolPwrss",tolPwrssvec[k],envir=environment(ff))
        t0 <- system.time(cur.opt <- try(optimx(val, ff,
                                                lower=c(0, rep(-Inf, 4L)),
                                                control=control,
                                                method="bobyqa")))
        assign("tolPwrss",tolPwrss0,envir=environment(ff))
        if (inherits(cur.opt,"try-error")) {
          d <- dbad
          d$result <- attr(cur.opt,"condition")$message
        } else {
          d <- d0
          d <- c(d,list(options=c(control,list(tolPwrss=tolPwrssvec[k])),
                        time=sum(t0[1:2]), ## user +system
                        parameters=cur.opt$par$par,           ## ??? why nested
                        deviance=c(cur.opt$fvalues$fvalues),  ## ??? why nested
                        KKT=c(cur.opt$KKT1,cur.opt$KKT2),
                        bad=NA,
                        result=NA_character_))
        }
        results[[ctr]] <- d
        save("results","sessinfo",file=fn)
      }
    }
  }
