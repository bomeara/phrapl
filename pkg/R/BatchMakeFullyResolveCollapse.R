BatchMakeFullyResolveCollapse<-function() {
source("../R/batchMS.R")
for (npop in 1:5) {
	for(maxK in 1:10) {
		print(paste("Now starting ",npop," populations and maxK=",maxK))
		start=proc.time()
		try(generateMigrationIndividualsFullyResolvedCollapseAllowNoMigration(popVector=rep(5,npop),maxK=maxK,verbose=TRUE,file=paste("migrationArray_FullyResolvedCollapse_npop",npop,"_maxK",maxK,".Rsave",sep="")))	
		print("That took")
		print(proc.time()-start)
	}
}
}