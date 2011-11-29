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
searchContinuousModelSpace(p=c(1,1,5), migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, nTrees=10000 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.tre",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=TRUE, debug=FALSE)
searchDiscreteModelSpace(migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, nTrees=10000 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.tre",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=TRUE, debug=FALSE,itnmax=2,method="BFGS",print.level=2)

#note: sometimes, with high migration rates (like 2.458274e+131) ms takes a really long time, sometimes needing to be quit
#this is something to deal with 