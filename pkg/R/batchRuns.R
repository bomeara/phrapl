library(parallel)
itnmaxVector<-c(2,25,100)
nTreesVector<-c(10000,1000000,100000000)
nTreesObservedVector<-c(3,10,100)
nrep<-10
doSingleRun<-function(itnmax,nTrees,nTreesObserved,rep) {
  newDir<-paste("itnmax",itnmax,"_nTrees",nTrees,"_nTreesObs",nTreesObserved,"_rep",rep,sep="")
  filename<-paste(newDir,".Rsave",sep="")
  system(paste("mkdir ",newDir))
  setwd(newDir)
  print(system("pwd"))
  system('cp ../*.pl .')
  source("../individualRunTestPerformance.R")
  individualRunTestPerformance(filename=filename,itnmax=itnmax, nTrees=nTrees, nTreesObserved=nTreesObserved)
}
singleParamDoSingleRun<-function(itnmax_nTrees_nTreesObserved_rep) {
  x<-as.integer(strsplit(itnmax_nTrees_rep,split="_")[[1]])
  doSingleRun(itnmax=x[1],nTrees=x[2],nTreesObserved=x[3],rep=x[4])
}
itnmax_nTrees_nTreesObservedVector<-c(outer(c(outer(paste(itnmaxVector,"_",sep=""),nTreesVector,"paste",sep="")),nTreesObservedVector,"paste",sep="_"))
itnmax_nTrees_nTreesObserved_repVector<-c(outer(itnmax_nTrees_nTreesObservedVector,sequence(nrep),"paste",sep="_"))

mc.cores<-max(1,detectCores()-2)
mclapply(itnmax_nTrees_nTreesObserved_repVector,FUN=singleParamDoSingleRun,mc.cores=mc.cores)
