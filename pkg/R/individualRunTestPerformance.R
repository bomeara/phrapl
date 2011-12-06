individualRunTestPerformance<-function(filename="testRunResult.Rsave",batchMSlocation="../batchMS.R",compareLocation="../comparecladespipe.pl",msLocation="/home/bomeara/bin/ms",popVector=c(3,4,4),maxK=2,itnmax=40,nTrees=1000000,migrationArray=NULL,migrationArrayMap=NULL, trueModelID=5,trueModelParams=c(0.9,0.2),nTreesObserved=3) {
  source(batchMSlocation)
  if(is.null(migrationArray)) {
    migrationArray<-generateMigrationIndividualsAllowNoMigration(popVector,maxK=maxK)
  }
  if(is.null(migrationArrayMap)) {
    migrationArrayMap<-generateMigrationArrayMap(migrationArray)
  }
  createAssignment(popVector)
  trueModel<-migrationArray[[trueModelID]]
  names(trueModelParams)<-msIndividualParameters(trueModel)
  msCallInfo<-createMSstringSpecific(popVector,trueModel,trueModelParams,nTrees=3)
  system(paste(msLocation=msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';' > observed.tre"),intern=FALSE)
  start.time<-Sys.time()
  result<-searchDiscreteModelSpaceOptim(migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, nTrees=10000 ,msLocation=msLocation,compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.tre",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=TRUE, debug=FALSE,itnmax=itnmax,method="BFGS",print.level=2)
  print(result)
  elapsed.time<-Sys.time()-start.time
  print(elapsed.time)
  save(result,elapsed.time,file=filename,compress=TRUE)
}