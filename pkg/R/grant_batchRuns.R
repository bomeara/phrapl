itnmaxVector<-c(2,25,100)
nTreesVector<-c(10000,1000000,100000000)
nTreesObservedVector<-c(3,10,100)
nrep<-10
dataDir<-"/Users/bomeara/Documents/MyDocuments/Active/phrapl/pkg/data/"
doSingleRun<-function(overallModel,itnmax,nTrees,nTreesObserved,rep,modelsToRemove,dataDir) {
  newDir<-paste("itnmax",itnmax,"_nTrees",nTrees,"_nTreesObs",nTreesObserved,"_rep",rep,sep="")
  filename<-paste(newDir,".Rsave",sep="")
  system(paste("mkdir ",newDir))
  setwd(newDir)
  print(system("pwd"))
  system('cp ../*.pl .')
  source("../grant_individualRunTestPerformance.R")
  popVector<-c(NA)
  trueModelParams<-c(NA)
  maxK<-2
  if (overallModel==1) { 
    #original plan
    #A four-population divergence model with some size change
  #For the phylogenetic model, lets do a pectinate tree with diversification occurring at 5N for the deepest split, then 3N and 1.5N. 
  #Theta = 5.0, but have it change of the branches to 2.5. We typically see values in this range for our IM runs in the salamanders.

    #modified to be three population divergence model
    #collapse 1 and 2, then 12 and 3, so tree is ((1,2),3)
    #     $collapseMatrix
    #      [,1] [,2]
    # [1,]    1    1
    # [2,]    1   NA
    # [3,]    0    1
    #Deepest split is 5N, most recent split is 3N
    #Have pop  1 with n0 multiplier of 5 after the split from 1
    #Have no migration between populations

    load(paste(dataDir,"migrationArray_npop3_maxK6.Rsave",sep=""))
    popVector<-rep(5,3)
  	maxK<-6 #need it at least 3 to get tree
    trueModelID<-2727
    modelsToRemove<-c() #models with migration
    for (i in sequence(length(migrationArray))) {
      if(sum(grepl("migration",msIndividualParameters(migrationArray[[i]])))!=0) {
        modelsToRemove<-append(modelsToRemove,i) 
      }
    }
    trueModelParams<-(3,5,5,1)
    names(trueModelparams)<-msIndividualParameters(migrationArray[[trueModelID]]
  }
  else if (overallModel==2) { #a divergence with gene flow model (2 pops)
  #For the divergence with gene flow, let theta-A be 5.0, theta1 be 3.0, and theta2 be 4.0. 
  #Divergence of 1.5N, and migration rates of 0.125.
  	popVector<-rep(5,2)
  	maxK<-4
  }
  else if (overallModel==3) { #a four population island model
  #For the n-island model, let theta1=6, theta2=6, theta3=4, and theta4=4. 
  #Symmetric migration matrix of 1-2 = 0.125, 1-3 = 0.250, 1-4 = 0.0625, 2-3 = 0, 2-4 = 0.125, and 3-4 = 0.250 
  #filter such that only migration models are possible -- no collapse
  	popVector<-rep(5,4)
  	maxK<-5
  
  }
  individualRunTestPerformance(filename=filename,itnmax=itnmax, nTrees=nTrees, nTreesObserved=nTreesObserved,modelsToRemove=modelsToRemove,maxK=maxK)
}
singleParamDoSingleRun<-function(itnmax_nTrees_nTreesObserved_rep,overallModel,modelsToRemove,dataDir) {
  x<-as.integer(strsplit(itnmax_nTrees_nTreesObserved_rep,split="_")[[1]])
  doSingleRun(overallModel=overallModel,itnmax=x[1],nTrees=x[2],nTreesObserved=x[3],rep=x[4],modelsToRemove=modelsToRemove,dataDir=dataDir)
}
itnmax_nTrees_nTreesObservedVector<-c(outer(c(outer(paste(itnmaxVector,"_",sep=""),nTreesVector,"paste",sep="")),nTreesObservedVector,"paste",sep="_"))
itnmax_nTrees_nTreesObserved_repVector<-c(outer(itnmax_nTrees_nTreesObservedVector,sequence(nrep),"paste",sep="_"))

mc.cores<-max(1,floor(detectCores()/3))
mclapply(itnmax_nTrees_nTreesObserved_repVector,FUN=singleParamDoSingleRun,mc.cores=mc.cores,dataDir=dataDir)
