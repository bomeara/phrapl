setwd("/Users/bomeara/Dropbox/CarstensScratch/analysis")
source("~/Documents/MyDocuments/Active/phrapl/pkg/R/batchMS.R")
itnmax=40
popVector<-c(3,4,4)
maxK<-2
migrationArray<-generateMigrationIndividualsAllowNoMigration(popVector,maxK=maxK)
migrationArrayMap<-generateMigrationArrayMap(migrationArray)
createAssignment(popVector)
trueModelID<-5
trueModel<-migrationArray[[trueModelID]]
trueModelParams<-c(0.9,0.2)
names(trueModelParams)<-msIndividualParameters(trueModel)
msCallInfo<-createMSstringSpecific(popVector,trueModel,trueModelParams,nTrees=3)
system(paste(msLocation="/usr/local/bin/ms",sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';' > observed.tre"),intern=FALSE)

returnAIC(c(log(1),log(1.5)),popVector,migrationIndividual=trueModel,nTrees=10000,compareLocation="comparecladespipe.pl", observed="observed.tre",print.results=TRUE,print.ms.string=TRUE,debug=TRUE)
searchContinuousModelSpaceOptim(p=c(1,1,5), migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, nTrees=10000 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.tre",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=TRUE, debug=FALSE,method="Nelder-Mead",itnmax=itnmax)
#searchContinuousModelSpace(p=c(1,1,5), migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, nTrees=10000 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.tre",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=TRUE, debug=FALSE,method="nlm",itnmax=2)
searchDiscreteModelSpaceOptim(migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, nTrees=10000 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.tre",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=TRUE, debug=FALSE,itnmax=itnmax,method="BFGS",print.level=2)
 
