## TO DO:
##   incorporate user, addr, etc. as appropriate args/options
##   alternative: (1) run from local machine, keep job open until finished
##                (2) run from SN
##   modularize: "submit single job"; "collect single job"; "clean up"
##   "check status" (sqjobs, parse output)
##   add options for wallclock, memory, etc.
## specialize SNjob
##  ability to restart (search for files, run relevant jobs)

## do scp transmission in one big block (faster than renegotiating for lots of small files)
## better clean-up / directory and file management


require(digest)  ## for job names

SNtarget <- "kraken"
SNaddr <- function() {
  paste(SNtarget,".sharcnet.ca",sep="")
}
SNuser <- "bolker"
SNfuser <- function() {
  paste(SNuser,"@",SNaddr(),sep="")
}
                                                   
scp <- function(from,to) {
  system(paste("scp",from,to))
}

sshcmd <- function(cmd,intern=FALSE) {
  system(paste("ssh",SNfuser(),cmd),intern=intern)
}

SNput <- function(fn,dir) {
  system(paste("ssh",SNaddr(),"mkdir -p",dir))
  scp(fn,paste(SNfuser(),":",dir,sep=""))
}

## use "" to use full path name for fn
SNget <- function(fn,dir) {
  scp(paste(SNfuser(),":",dir,"/",fn,sep=""),getwd())
}

sqcmd <- function(cmd,intern=FALSE,verbose=FALSE) {
  res <- character(length(cmd))
  for (i in seq_along(cmd)) {
    scmd <- paste("'bash -l -c \"sq",cmd[i],"\"\'",sep="")
    if (verbose) cat(scmd,"\n")
    res[i] <- system(paste("ssh",SNaddr(),scmd),
                     intern=intern,wait=TRUE)
  }
  res
}

sqjobs <- function() sqcmd("jobs")

SNjob <- function(batchname,otherfiles,wtime="2h",
                  mem="2G",SNwd="/work/bolker/projects/pines",
                  verbose=FALSE) {
  b <- c(paste("setwd('",SNwd,"')",sep=""),readLines(paste(batchname,".R",sep="")))
  SNbatchfile <- paste(batchname,"_SN.R",sep="")
  writeLines(b,con=SNbatchfile)
  SNput(SNbatchfile,"~")
  sapply(otherfiles,SNput,SNwd)
  ss <- sqcmd(paste("sub -j",batchname,"--mpp",mem,
        "-r",wtime,
        "-o",paste(batchname,"out",sep="."),
        "R CMD BATCH",SNbatchfile),intern=TRUE,verbose=verbose)
  SNjobid <- gsub("submitted as jobid ","",ss)
  SNjobid
} 

## retrieve files
SNgetpat <- function(pat,wd) {
  files <- sshcmd(paste("ls ",wd,"/",pat,sep=""),intern=TRUE)
  sapply(files,SNget,"")
}

SNlapply <- function(X,FUN,clean=TRUE,
                     Robjs=NULL,Rfiles=NULL,
                     mem="2G",wtime="30m") {
  jobid <- digest(c(Sys.time(),SNuser),algo="crc32")
  SNwd <- file.path("/work",SNuser,paste("Rjobs",jobid,sep="_"))
  sshcmd(paste("mkdir -p",SNwd))
  ## write individual-chunk data
  chunkids <- character(length(X))
  jvec <- seq_along(X)
  chunkids <- sprintf("%s_%03d",jobid,jvec)
  tt <- paste(chunkids,"RData",sep=".")
  mapply(function(x,f) {
    save("x",file=f) ## chunk-specific data saved as 'x'
  },X,tt)
  sapply(tt,SNput,SNwd)
  rm(X)               ## clean-up in case we want to save the whole image somehow
  ttlist_out <- paste(chunkids,"_out",".RData",sep="")
  ttlist_R <- paste(chunkids,".R",sep="")
  do.call(save,c(list("FUN"),as.list(Robjs),list(file="image.RData")))  ## save FUN and required objects
  sapply(c(as.list(Rfiles),"image.RData"),SNput,SNwd)
  ppaste <- function(...) paste(...,sep="")
  Rcmds <- paste(ppaste("setwd('",SNwd,"');"),
                 "load('image.RData');",
                 ppaste("load ('",tt,"');"),
                 "y <- FUN(x);",
                 ppaste("save('y',file='",ttlist_out,"')"))
  invisible(mapply(writeLines,as.list(Rcmds),as.list(ttlist_R)))
  sapply(ttlist_R,SNput,"~")
  spawncmds <- paste("sub -j",chunkids,
                     "--mpp",mem,
                     "-r",wtime,
                     "-o",paste(chunkids,"out",sep="."),
                     "R CMD BATCH",ttlist_R)
  submitstr <- sqcmd(spawncmds,intern=TRUE)
  SNjobids <- gsub("submitted as jobid ","",submitstr)
  ## run --waitfor job
  wjobs <- paste("-w",paste(SNjobids,collapse=","))
  collectname <- paste("collect_",jobid,sep="")
  collectfile <- paste(collectname,".R",sep="")
  alldatafile <- paste(jobid,"_ALL.RData",sep="")
  ## read in files (note return variable is called 'y' in every instance)
  writeLines(c(ppaste("setwd('",SNwd,"')"),
               paste("ttlist_out <- list.files(pattern='",jobid,"_[0-9]+_out.RData')",sep=""),
               "Y <- lapply(ttlist_out,function(x) { load(x); y })",
               paste("save('Y',file='",alldatafile,"')",sep="")),
             con=collectfile)
  SNput(collectfile,"~")
  sqcmd(paste("sub -j",collectname,
              wjobs,
              "--mpp 2G",
              "-r 10m",
              "-o",paste(collectname,"out",sep="."),
              "R CMD BATCH",collectfile))
  ## now need to wait until it's finished ...
  get_results <- function() {
    SNget(alldatafile,SNwd)
    load(alldatafile)
    Y
  }
  ## clean up?
  ## if (clean) unlink(c(ttlist,ttlist_out))
  unlink(list.files(pattern=paste(jobid,"_[0-9]+",sep="")))
  ## cleans all but "collect" and ALL.RData
  list(get_results)
}

  if (FALSE) {
## tests
  z <- 1:5
  fun <- function(x) x^2
  debug(SNlapply)
  zz <- SNlapply(as.list(z),fun)
}
