setwd("/Users/bomeara/Dropbox/InProgress/Carstens/Basic")
popVector<-rep(7,3)
migrationArray<-generateMigrationIndividuals(popVector)
createAssignment(popVector)
allValues<-c(5,.2,.02,.0001,.001,.01)
names(allValues)<-c("collapse_1","n0multiplier_1","n0multiplier_2","migration_1","migration_2","migration_3")
print(allValues)
trueModel<-57
nTrees=10000
#trueModel<-1
numModels<-length(migrationArray)
saveMS(popVector,migrationArray[[trueModel]],allValues,nTrees=25,file="obs.tre")
results<-data.frame()
for (model in 1:numModels) {
  saveMS(popVector,migrationArray[[model]],allValues,nTrees=nTrees,file="sim.tre")
  output<-system("perl compareclades.pl -aassign.txt -oobs.tre -ssim.tre",intern=TRUE)
  lnL<-convertOutputVectorToLikelihood(output,nTrees)
  AIC<-2*length(msIndividualParameters(migrationArray[[model]])) - 2*lnL
  print(paste(model,lnL,AIC,paste(msIndividualParameters(migrationArray[[model]]),sep=" ",collapse=" ")))
  if (model==1) {
   results<-data.frame(model,lnL,AIC, paste(msIndividualParameters(migrationArray[[model]]),sep=" ",collapse=" "))
  }
  else {
   results<-rbind(results,data.frame(model,lnL,AIC, paste(msIndividualParameters(migrationArray[[model]]),sep=" ",collapse=" ")))   
  }
}
results<-cbind(trueModel=(results$model==trueModel),deltaAIC=results$AIC-min(results$AIC),results)
print(results)