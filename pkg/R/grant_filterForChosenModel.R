
#for model1
load(paste(dataDir,"migrationArray_npop3_maxK6.Rsave",sep=""))
    

goodcount<-0
modelsToRemove<-c() #models with migration
for (i in sequence(length(migrationArray))) {
  if(sum(grepl("migration",msIndividualParameters(migrationArray[[i]])))!=0) {
    modelsToRemove<-append(modelsToRemove,i) 
  }
    if(sum(grepl("collapse",msIndividualParameters(migrationArray[[i]])))==2)  {
    
    if(sum(grepl("n0multiplier",msIndividualParameters(migrationArray[[i]])))==2)  {
          if(sum(grepl("migration",msIndividualParameters(migrationArray[[i]])))==0)  {

    print("new one")
    print(i)
    print(msIndividualParameters(migrationArray[[i]]))
    print(sum(grepl("migration",msIndividualParameters(migrationArray[[i]]))))
    print(migrationArray[[i]])
    goodcount<-goodcount+1
          }
    }
  }
}
print(goodcount)


