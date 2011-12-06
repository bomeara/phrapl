itnmaxVector<-c(2,25,100)
nTreesVector<-c(10000,1000000,100000000)
nrep<-6
doSingleRun<-function(itnmax,nTrees,rep) {
  newDir<-paste("itnmax",itnmax,"_nTrees",nTrees,"_rep",rep,sep="")
  filename<-paste(newDir,".Rsave",sep="")
  system(paste("mkdir ",newDir))
  setwd(newDir)
  load("../individualRunTestPerformance.R")
  individualRunTestPerformance(filename=filename,itnmax=itnmax, nTrees=nTrees)
}
singleParamDoSingleRun<-function(itnmax_nTrees_rep) {
  x<-as.integer(strsplit(itnmax_nTrees_rep,split="_")[[1]])
  doSingleRun(itnmax=x[1],nTrees=x[2],rep=x[3])
}
itnmax_nTrees_repVector<-c(outer(c(outer(paste(itnmaxVector,"_",sep=""),nTreesVector,"paste",sep="")),sequence(nrep),"paste",sep="_"))

mc.cores<-detectCores()
mclapply(itnmax_nTrees_repVector,FUN=singleParamDoSingleRun,mc.cores=mc.cores)