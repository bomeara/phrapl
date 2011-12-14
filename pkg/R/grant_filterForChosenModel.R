
#for model3
load(paste(dataDir,"migrationArray_CollapseForbidden_npop4_maxK6.Rsave",sep=""))
migrationArray<-returnSymmetricMigrationOnly(migrationArray)


goodcount<-0
for (i in sequence(length(migrationArray))) {
    print(cat("n0",sum(grepl("n0multiplier",msIndividualParameters(migrationArray[[i]]))),"mig", sum(grepl("migration",msIndividualParameters(migrationArray[[i]])))))
    if(sum(grepl("n0multiplier",msIndividualParameters(migrationArray[[i]])))==3)  {
          if(sum(grepl("migration",msIndividualParameters(migrationArray[[i]])))==3)  {

    print("new one")
    print(i)
    #print(msIndividualParameters(migrationArray[[i]]))
    print(sum(grepl("migration",msIndividualParameters(migrationArray[[i]]))))
    print(migrationArray[[i]])
    goodcount<-goodcount+1
          }
    }
  }

print(goodcount)


