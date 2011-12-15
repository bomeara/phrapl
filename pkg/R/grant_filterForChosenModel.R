
#for model3
load(paste(dataDir,"migrationArray_CollapseForbidden_npop4_maxK6.Rsave",sep=""))
trueModel<-migrationArray[[1]]
trueModel$n0multiplierMap<-matrix(c(1,1,2,2),ncol=1)
trueModel$migrationArray[, , 1]<-matrix(c(NA,1,2,3,NA,NA,0,1,NA,NA,NA,2,NA,NA,NA,NA),ncol=4,byrow=TRUE)
migrationArray[[length(migrationArray)+1]]<-trueModel
migrationArray<-returnSymmetricMigrationOnly(migrationArray)


goodcount<-0
for (i in sequence(length(migrationArray))) {
    print(cat("n0",sum(grepl("n0multiplier",msIndividualParameters(migrationArray[[i]]))),"mig", sum(grepl("migration",msIndividualParameters(migrationArray[[i]])))))
    focalArray<-migrationArray[[i]]$migrationArray[ , , 1]
    print(focalArray[upper.tri(focalArray)])
    if(max(abs(focalArray[upper.tri(focalArray)]-c(1, 2, 3, 0, 1, 2)))==0) {
    
  #  if(sum(grepl("n0multiplier",msIndividualParameters(migrationArray[[i]])))==3)  {
  #        if(sum(grepl("migration",msIndividualParameters(migrationArray[[i]])))==3)  {

    print("new one")
    print(i)
    #print(msIndividualParameters(migrationArray[[i]]))
    print(sum(grepl("migration",msIndividualParameters(migrationArray[[i]]))))
    print(migrationArray[[i]])
    goodcount<-goodcount+1
          }
   # }
  }

print(goodcount)


