donnes <- read.csv("donnesR_min.csv")
library(lme4Eigen)
sessionInfo()$otherPkgs$lme4Eigen$Version

form <- coupe ~ diffdhp + marchand:(diffdhp + diffdhp2 + MSCR_ID) +
  essence + (1|IDENT)

optimizer <- c("bobyqa")
nAGQ <- c(1:5,8,12,20,25)

## loop over tolPwrss, rhobeg, rhoend?

d0 <- list(problem="delisle_forest")

tmpf <- function(optimizer,nAGQ) {
  cur.opt <- glmer(form,donnes,binomial,
                   optimizer=optimizer,nAGQ=nAGQ)
  if (inherits(cur.opt,"try-error")) {
    d <- dbad
    d$result <- attr(cur.opt,"condition")$message
  } else {
    d <- d0
    d <- c(d,list(method=optimizer,
                  options=NA,
                  nAGQ=nAGQ,
                  time=sum(t0[1:2]), ## user +system
                  parameters=allcoef(cur.opt),
                  deviance=deviance(cur.opt),
                  KKT=c(NA,NA),
                  bad=NA,
                  result=NA_character_))
  }
  d
}

                   
tmpf("bobyqa",1)
tmpf("NelderMead",1)

xapply(
t1 <- system.time(MSCR.mixed01 <- glmer(form, donnes, binomial, nAGQ=1L))
allcoef(MSCR.mixed01)
fixef(MSCR.mixed01)

MSCR.mixed <- glmer(coupe ~ diffdhp + marchand:diffdhp + marchand:diffdhp2 + essence + marchand:MSCR_ID + (1|IDENT), family = binomial)
summary(MSCR.mixed)
