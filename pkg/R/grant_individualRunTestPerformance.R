grant_individualRunTestPerformance<-function(filename="testRunResult.Rsave",batchMSlocation="../batchMS.R",compareLocation="comparecladespipe.pl",msLocation="/usr/local/bin/ms",popVector=c(3,4,4),maxK=2,itnmax=40,nTrees=1000000,migrationArray=NULL,migrationArrayMap=NULL, trueModelID=5,trueModelParams=c(0.9,0.2),nTreesObserved=3,modelsToRemove=NULL) {
  source(batchMSlocation)
  print("creating migrationArray")
  if(is.null(migrationArray)) {
    migrationArray<-generateMigrationIndividualsAllowNoMigration(popVector,maxK=maxK)
  }
  print("getting true model")
  trueModel<-migrationArray[[trueModelID]] #NOTE: We use a tree model from the whole set before tossing modelsToRemove
  #This lets us use a model for generation we don't have in the analysis set

  if(!is.null(modelsToRemove)) {
  	migrationArray<-migrationArray[-modelsToRemove]
  }
  if(is.null(migrationArrayMap) || !is.null(modelsToRemove)) {
    migrationArrayMap<-generateMigrationArrayMap(migrationArray)
  }
  print(paste("number of models in migrationArray = ",length(migrationArray)))
  createAssignment(popVector)
  names(trueModelParams)<-msIndividualParameters(trueModel)
  msCallInfo<-createMSstringSpecific(popVector,trueModel,trueModelParams,nTrees=nTreesObserved)
  print(paste("ms call is ",msCallInfo))
  system(paste(msLocation=msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';' > observed.tre"),intern=FALSE)
  start.time<-proc.time()
  result<-searchDiscreteModelSpaceOptim(migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, nTrees=10000 ,msLocation=msLocation,compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.tre",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=TRUE, debug=TRUE,itnmax=itnmax,method="BFGS",print.level=2)
  print(result)
  elapsed.time<-proc.time()-start.time
  print(elapsed.time)
  recoveredModelID<-returnModel(result$par,migrationArrayMap)
  summaryVec<-c(trueModelID,recoveredModelID,itnmax,nTrees,nTreesObserved,elapsed.time)
  names(summaryVec)<-c("trueModelID","recoveredModelID","itnmax","nTrees","nTreesObserved","seconds")
  print(summaryVec)
  cat(names(summaryVec),file="summary.txt")
  cat("\n",file="summary.txt",append=TRUE)
  
  cat(summaryVec,file="summary.txt",append=TRUE)
  cat("\n",file="summary.txt",append=TRUE)
  cat("\nelapsed.time",file="summary.txt",append=TRUE)
  cat(elapsed.time,file="summary.txt",append=TRUE)
  cat("\ncollapse matrix",file="summary.txt",append=TRUE)
  cat("\ntrue  ",file="summary.txt",append=TRUE)
  cat(migrationArray[[trueModelID]]$collapseMatrix,file="summary.txt",append=TRUE)
  cat("\ninfer ",file="summary.txt",append=TRUE)
  cat(migrationArray[[recoveredModelID]]$collapseMatrix,file="summary.txt",append=TRUE)

  cat("\n\nn0multiplierMap",file="summary.txt",append=TRUE)
  cat("\ntrue  ",file="summary.txt",append=TRUE)
  cat(migrationArray[[trueModelID]]$n0multiplierMap,file="summary.txt",append=TRUE)
  cat("\ninfer ",file="summary.txt",append=TRUE)
  cat(migrationArray[[recoveredModelID]]$n0multiplierMap,file="summary.txt",append=TRUE)

  cat("\n\nmigrationArray",file="summary.txt",append=TRUE)
  cat("\ntrue  ",file="summary.txt",append=TRUE)
  cat(migrationArray[[trueModelID]]$migrationArray,file="summary.txt",append=TRUE)
  cat("\ninfer ",file="summary.txt",append=TRUE)
  cat(migrationArray[[recoveredModelID]]$migrationArray,file="summary.txt",append=TRUE)

  
  save(list=ls(),file=filename,compress=TRUE)
  
}
