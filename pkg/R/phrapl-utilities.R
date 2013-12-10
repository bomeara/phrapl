#library(partitions)
#library(ape)
#library(rgenoud)
#library(optimx)
#library(parallel)
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ncores", "popVector", "maxK","migrationArray", "migrationArrayMap"))

ColMax<-function(x,na.rm=TRUE) {
	maxVal=rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		maxVal[i]<-max(x[,i],na.rm=na.rm)
	}
	return(maxVal)
}

ColMin<-function(x,na.rm=TRUE) {
	minVal=rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		minVal[i]<-min(x[,i],na.rm=na.rm)
	}
	return(minVal)
}

RowMax<-function(x,na.rm=TRUE) {
  maxVal=rep(NA,dim(x)[1])
	for (i in 1:dim(x)[1]) {
		maxVal[i]<-max(x[i,],na.rm=na.rm)
	}
	return(maxVal)
}

RowMin<-function(x,na.rm=TRUE) {
	minVal=rep(NA,dim(x)[1])
	for (i in 1:dim(x)[1]) {
		minVal[i]<-min(x[i,],na.rm=na.rm)
	}
	return(minVal)
}

MatrixContainsAllValues<-function(x,minVal=1) {
	maxVal<-max(x,na.rm=TRUE)
	hasAllValues<-TRUE
	for(i in c(minVal,maxVal)) {
		if(sum(x==i,na.rm=TRUE)==0) {
			hasAllValues<-FALSE
			return(hasAllValues)
		}
	}
	return(hasAllValues)
}

ColCountIf<-function(x,val) {
	countVec<-rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		countVec[i]=length(which(x[,i]==val))
	}
	return(countVec)
}

ColCountLast<-function(x,val) {
	return(length(which(x[,dim(x)[2]]==val)))
}



Popinterval<-function(collapseMatrix,complete=FALSE) {
#collapse matrix has populations as rows and each generation going back in time as columns
	pInt<-list(collapseMatrix=collapseMatrix,complete=complete)
	class(pInt)="popinterval"
	if (max(collapseMatrix,na.rm=TRUE)==0) {
		pInt$complete=TRUE #if the model is one of no collapse, nothing further to do
	}
	return(pInt)
}

N0multiplierindividual<-function(collapseMatrix,complete=FALSE,n0multiplierMap) {
	tI<-list(collapseMatrix=collapseMatrix,complete=complete,n0multiplierMap=n0multiplierMap)
	class(tI)<-"n0multiplierindividual"
	return(tI)
}

Migrationindividual<-function(collapseMatrix,complete=FALSE,n0multiplierMap,migrationArray) {
	tI<-list(collapseMatrix=collapseMatrix,complete=complete,n0multiplierMap=n0multiplierMap,migrationArray=migrationArray)
	class(tI)<-"migrationindividual"
	return(tI)
}

Expansionmultiplierindividual<-function(collapseMatrix, complete, expansionmultiplierMap) {
	tI<-list()  #FIX THIS!!!!
	stop("This code is not written yet")
	return(tI)
}



CreateMSstringSpecific<-function(popVector,migrationIndividual,parameterVector,nTrees=1) {
	collapseMatrix<-migrationIndividual$collapseMatrix
	complete<-migrationIndividual$complete
	n0multiplierMap<-migrationIndividual$n0multiplierMap
	migrationArray<-migrationIndividual$migrationArray
	nPop<-length(popVector)
	numSteps<-dim(collapseMatrix)[2]
	n0multiplierParameters<-parameterVector[grep("n0multiplier",names(parameterVector))]
	migrationParameters<-parameterVector[grep("migration",names(parameterVector))]
	collapseParameters<-parameterVector[grep("collapse",names(parameterVector))]
	
	if (1>2 ) {
# stop(paste("Incorrect parameterVector: you passed ",length(parameterVector),"entries but it needs",length(parameterList),"entries:",paste(parameterList,sep=" ",collapse="")))
	}
	else {
		msString<-paste("-T -I",nPop,paste(popVector,sep=" ", collapse=" "),sep=" ")
		
#do values at present
		initialN0multipler<-""
		for (i in 1:nPop) {
			initialN0multipler<-paste(initialN0multipler,"-n",i,sprintf("%f",n0multiplierParameters[n0multiplierMap[i,1] ]),sep=" ")
		}
		
		initialMigration<-"-ma"
		for (i in 1:nPop) {
			for (j in 1:nPop) {
				if (i==j) {
					initialMigration<-paste(initialMigration, "x", sep=" ") 
				}
				else {
					if (migrationArray[i,j,1] > 0) {
						initialMigration<-paste(initialMigration, sprintf("%f",migrationParameters[migrationArray[i,j,1] ]),sep=" ")
					}
					else {
						initialMigration<-paste(initialMigration, sprintf("%f",0),sep=" ")
					}
				}
			}
		}
		msString<-paste(msString,initialN0multipler,initialMigration,sep=" ")
		
#is there a collapse in the first gen?
		if (max(collapseMatrix[,1],na.rm=TRUE)>0) { #remember that everything goes to the lowest index population with state ==1
			initialCollapse<-""
			popsToCollapse<-which(collapseMatrix[,1]==1)
			for (i in 2:length(popsToCollapse)) {
				initialCollapse<-paste(initialCollapse, "-ej", sprintf("%f",collapseParameters[1]), popsToCollapse[i], popsToCollapse[1]  ,sep=" ")       
			}
			msString<-paste(msString,initialCollapse,sep=" ")
		}
		
#now go back in time
		if (numSteps>1) {
			for (generation in 2:numSteps) {
				collapseTime<-collapseParameters[generation-1]
				n0multiplierPositions<-which(!is.na(n0multiplierMap[,generation]))
				for (n0multiplierPosIndex in 1:length(n0multiplierPositions)) {
					n0multiplierPos<-n0multiplierPositions[n0multiplierPosIndex]
					msString<-paste(msString,"-en",sprintf("%f",collapseTime),sprintf("%f",n0multiplierPos),sprintf("%f",n0multiplierParameters[n0multiplierMap[n0multiplierPos,generation] ]),sep=" ")
				}
				
				for (fromPos in 1:nPop) {
					for (toPos in 1:nPop) {
						migrationArrayValue<-migrationArray[fromPos,toPos,generation]
						if (!is.na(migrationArrayValue) ) {
							if (migrationArrayValue>0) {
								msString<-paste(msString,"-em",sprintf("%f",collapseTime),fromPos,toPos,sprintf("%f",migrationParameters[migrationArrayValue]),sep=" ")     
							}
							else {
								msString<-paste(msString,"-em",sprintf("%f",collapseTime),fromPos,toPos,sprintf("%f",0),sep=" ")     
							}
						}
					}
				}
				
				
#is there a collapse in this gen?
				if (max(collapseMatrix[,generation],na.rm=TRUE)>0) { #remember that everything goes to the lowest index population with state ==1
					initialCollapse<-""
					popsToCollapse<-which(collapseMatrix[,generation]==1)
					for (i in 2:length(popsToCollapse)) {
						initialCollapse<-paste(initialCollapse, "-ej", sprintf("%f",collapseParameters[generation]), popsToCollapse[i], popsToCollapse[1]  ,sep=" ")       
					}
					msString<-paste(msString,initialCollapse,sep=" ")
				}
				
			}
		}
		nsam<-sum(popVector)
		nreps<-nTrees
		opts<-msString
		returnobject<-list(nsam=nsam,nreps=nreps,opts=opts)
		return(returnobject)
	}
}

CreateMSstringGeneric<-function(popVector,migrationIndividual,parameterVector,nTrees=1) {
	return(CreateMSstringSpecific(popVector,migrationIndividual,MsIndividualParameters(migrationIndividual),nTrees)) 
}

ReturnIncompletes<-function(popIntervalsList) {
	incompleteElements<-c()
	for (i in 1:length(popIntervalsList)) {
		if (!popIntervalsList[[i]]$complete) {
			incompleteElements<-append(incompleteElements,i)
		}
	}
	return(incompleteElements)
}

UpdateCompletes<-function(popIntervalsList) {
	for (i in 1:length(popIntervalsList)) {
		if (ColCountLast(popIntervalsList[[i]]$collapseMatrix,1)==0) {
			popIntervalsList[[i]]$complete<-TRUE #end up with having no more collapses
		}
		else if (ColCountLast(popIntervalsList[[i]]$collapseMatrix,0)==0) {
			popIntervalsList[[i]]$complete<-TRUE #end up as one population
		}
	}
	return(popIntervalsList)	
}

CompleteIntervals<-function(popIntervalsList) {
	while(length(ReturnIncompletes(popIntervalsList))>0) {
		incompleteElements <- ReturnIncompletes(popIntervalsList)
		newPopIntervalsList <- popIntervalsList[-1*incompleteElements] #has only complete elements
		for (index in 1:length(incompleteElements)) {
			currentParent=popIntervalsList[[incompleteElements[index] ]]
			lastGen<-currentParent$collapseMatrix[,dim(currentParent$collapseMatrix)[2] ]
#all those in the lastGen in state 1 were merged into one population, identified by the one with the lowest id
			survPop=length(which(lastGen==0))
			if (length(which(lastGen==1))>0) {
				survPop=survPop+1
			}
			else {
				stop("this population was complete")
			}
			rawIntervals<-c()
			if (survPop>1) {
				rawIntervals<-blockparts(c(1:survPop),survPop,include.fewer=TRUE)
				rawIntervals<-rawIntervals[,ColMax(rawIntervals)<=1] #we are okay with having all zeros: no population collapse
				rawIntervals<-rawIntervals[,ColCountIf(rawIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
			}
			else {
				rawIntervals=matrix(c(0,1),nrow=1)
			}
			rowMapping<-c(min(which(lastGen==1)),which(lastGen==0))
			rawIntervalsRescaled<-matrix(NA,ncol=dim(rawIntervals)[2],nrow=length(lastGen))
			for (j in 1:length(rowMapping)) {
				rawIntervalsRescaled[rowMapping[j], ]<-rawIntervals[j,]
			}
			for (k in 1:dim(rawIntervalsRescaled)[2]) {
				newPopIntervalsList[[length(newPopIntervalsList)+1]]<-Popinterval(cbind(currentParent$collapseMatrix,matrix(rawIntervalsRescaled[,k],ncol=1)))
			}
		}
		popIntervalsList<-newPopIntervalsList
		popIntervalsList<-UpdateCompletes(popIntervalsList)
	}
	return(popIntervalsList)
}


GenerateIntervals<-function(popVector,maxK) {
#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
#maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
	nPop <- length(popVector)
	firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	firstIntervals<-firstIntervals[,ColMax(firstIntervals)<=1] #we are okay with having all zeros: no population collapse
	firstIntervals<-firstIntervals[,ColCountIf(firstIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
	popIntervalsList<-list()
	for (i in 1:dim(firstIntervals)[2]) {
		popIntervalsList[[i]]<-Popinterval(as.matrix(firstIntervals[,i],ncol=1))
	}
	popIntervalsList<-CompleteIntervals(UpdateCompletes(popIntervalsList))
	return(popIntervalsList)
}

GenerateIntervalsNoCollapse<-function(popVector,maxK) {
#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
#maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
  nPop <- length(popVector)
	firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	firstIntervals<-firstIntervals[,ColMax(firstIntervals)==0] #we only want ones with no collapse (like an island model)
	popIntervalsList<-list()
	popIntervalsList[[1]]<-Popinterval(as.matrix(firstIntervals,ncol=1))
	popIntervalsList<-CompleteIntervals(UpdateCompletes(popIntervalsList))
	return(popIntervalsList)
}

GenerateIntervalsFullyResolvedCollapse<-function(popVector,maxK) {
#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
#maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
  nPop <- length(popVector)
	firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	firstIntervals<-firstIntervals[,ColMax(firstIntervals)<=1] #we are okay with having all zeros: no population collapse
	firstIntervals<-firstIntervals[,ColCountIf(firstIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
	popIntervalsList<-list()
	for (i in 1:dim(firstIntervals)[2]) {
		popIntervalsList[[i]]<-Popinterval(as.matrix(firstIntervals[,i],ncol=1))
	}
	popIntervalsList<-CompleteIntervals(UpdateCompletes(popIntervalsList))
  intervalsToDelete<-c()
  for (i in sequence(length(popIntervalsList))) {
    focalCollapseMatrix<-popIntervalsList[[i]]$collapseMatrix
    if ((min(ColMax(focalCollapseMatrix))==0) || (dim(focalCollapseMatrix)[1]!=(1+dim(focalCollapseMatrix)[2]))) {
       intervalsToDelete<-append(intervalsToDelete,i)
    }
  }
  popIntervalsList<-popIntervalsList[-intervalsToDelete]
	return(popIntervalsList)
}


#the basic idea here is that at each population in each time interval there is a n0multiplier. These can all be set to the same value, allowed to vary, or assigned in clumps (i.e., pops 1, 3, and 6 have the same n0multiplier value)
#this generates all such mappings, subject to staying within the maximum number of free parameters
GenerateN0multiplierIndividuals<-function(popVector,popIntervalsList=GenerateIntervals(popVector,maxK=SetMaxK(popVector)),maxK=SetMaxK(popVector)) {
	n0multiplierIndividualsList<-list()
	for (i in 1:length(popIntervalsList)) {
		n0multiplierMapTemplate<-1+0*popIntervalsList[[i]]$collapseMatrix  #will have all the populations, all with either NA or 1
		numLineages=sum(n0multiplierMapTemplate,na.rm=TRUE)
		possibleMappings<-compositions(numLineages)
		for (mappingIndex in 1:dim(possibleMappings)[2]) {
			thisMapping<-possibleMappings[,mappingIndex]
			if ((length(which(thisMapping>0))+KPopInterval(popIntervalsList[[i]]) )<=maxK) { #only do it on those mappings that have not too many free parameters
				n0multiplierMap<-n0multiplierMapTemplate
				whichPositions <- which(n0multiplierMap==1)
				for (positionIndex in 1:length(whichPositions)) {
					position=whichPositions[positionIndex]
					paramPosition<-which(thisMapping>0)[1]
					n0multiplierMap[position]=paramPosition #the position of the first parameter
					thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we have used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
				}
				n0multiplierIndividualsList[[length(n0multiplierIndividualsList)+1]]<-N0multiplierindividual(popIntervalsList[[i]]$collapseMatrix, popIntervalsList[[i]]$complete, n0multiplierMap)
			}
		}
	}
	return(n0multiplierIndividualsList)
}

#now we will generate all possible assignments of pairwise migration. Again, we want to keep the total number of free parameters (times, n0multipliers, migration rates) under our chosen max
#allow a model where migrations change anywhere along branch, or only at coalescent nodes? The problem with the latter is that you cannott fit some reasonable models: i.e., two populations persisting through time. Problem with the former is parameter space
GenerateMigrationIndividuals<-function(popVector,n0multiplierIndividualsList=GenerateN0multiplierIndividuals(popVector,popIntervalsList=GenerateIntervals(popVector,maxK=SetMaxK(popVector)),maxK=SetMaxK(popVector)), maxK=SetMaxK(popVector), verbose=FALSE) {
	migrationIndividualsList<-list()
	for (i in sequence(n0multiplierIndividualsList)) {
		if (verbose==TRUE) {
			print(paste("doing ",i,"/",length(n0multiplierIndividualsList)))
		}
		collapseMatrix<-n0multiplierIndividualsList[[i]]$collapseMatrix
		n0multiplierMap<-n0multiplierIndividualsList[[i]]$n0multiplierMap
		numFinalPops<-dim(collapseMatrix)[1]
		numSteps<-dim(collapseMatrix)[2]
		if ((KCollapseMatrix(collapseMatrix) + KN0multiplierMap(n0multiplierMap) )<=maxK) {
			migrationTemplate<-array(data=NA,dim=c(numFinalPops,numFinalPops,numSteps),dimnames=c("from","to","generation"))
			for (interval in 1:numSteps) {
				for (fromPop in 1:numFinalPops) {
					for (toPop in 1:numFinalPops) {
						if (fromPop!=toPop) {
							if (!is.na(collapseMatrix[fromPop,interval]) && !is.na(collapseMatrix[toPop,interval])) {
								migrationTemplate[fromPop,toPop,interval]<-1
							}
						}
					}
				}
			}
			numPairs=sum(migrationTemplate,na.rm=TRUE)
		#	print(paste("starting compositions for ",numPairs))
		#	possibleMappings<-compositions(numPairs)
		#	print(paste("length of that was ",dim(possibleMappings)[2]))
		#	maxMappingK<-dim(possibleMappings)[1]
		#	actualMappingK <- maxMappingK - colCountIf(possibleMappings,0)
		#	netK<-actualMappingK + kCollapseMatrix(collapseMatrix) + kN0multiplierMap(n0multiplierMap)
			mappingMaxKAllowed <- maxK - ( KCollapseMatrix(collapseMatrix) + KN0multiplierMap(n0multiplierMap) )
		#	possibleMappings<-possibleMappings[,which(netK<=maxK)]
			
			evaluateMapping<-TRUE
			if (mappingMaxKAllowed<=0) {
				evaluateMapping<-FALSE #since there's no way to do this and not have too many parameters
			}
			thisMapping<-firstcomposition(numPairs)
			#for (mappingIndex in sequence(dim(possibleMappings)[2])) {
			while(evaluateMapping) {
				thisMapping<-c(thisMapping,rep(0,numPairs-length(thisMapping))) #adds trailing zeros
				thisMappingOrig<-thisMapping
				if ( length(which(thisMapping>0)) <= mappingMaxKAllowed) { #only do it on those mappings that have not too many free parameters
					migrationArray<-migrationTemplate
					whichPositions <- which(migrationArray==1)
					for (positionIndex in 1:length(whichPositions)) {
						position=whichPositions[positionIndex]
						paramPosition<-which(thisMapping>0)[1]
						migrationArray[position]=paramPosition #the position of the first parameter
						thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we have used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
					}
					migrationIndividualsList[[length(migrationIndividualsList)+1]]<-Migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArray)
					#print(paste("Just created migration individual ",length(migrationIndividualsList)))
				}
				if (!islastcomposition(thisMappingOrig,restricted=FALSE)) {
					thisMapping<-nextcomposition(thisMappingOrig,restricted=FALSE)
				}
				else {
					evaluateMapping<-FALSE
				}
			}		
		}
	}
	return(migrationIndividualsList)
}

#this will create a set of models with no populations merging. However, not all populations need to have migration to every other population
GenerateMigrationIndividualsNoCollapseAllowNoMigration<-function(popVector,n0multiplierIndividualsList=GenerateN0multiplierIndividuals(popVector,popIntervalsList=GenerateIntervalsNoCollapse(popVector,maxK=SetMaxK(popVector)),maxK=SetMaxK(popVector)), maxK=SetMaxK(popVector), verbose=FALSE,file=NULL) {
  migrationArray<-GenerateMigrationIndividualsAllowNoMigration(popVector,n0multiplierIndividualsList, maxK, verbose=verbose, file=NULL)
  migrationIndividualsToKill<-c()
  for (i in sequence(length(migrationArray))) {
     if (RowMax(migrationArray[[i]]$migrationArray[, , 1])==0 && RowMax(migrationArray[[i]]$migrationArray[, , 1])==0) {
        migrationIndividualsToKill<-append( migrationIndividualsToKill,i)
     }
  }
  migrationArray<-migrationArray[-1*migrationIndividualsToKill]
  if (!is.null(file)) {
    save(migrationArray,maxK,popVector,file=file,compress=TRUE)
  }
  return(migrationArray)
}

#this will create a set of models where populations are linked only by a bifurcating tree
GenerateMigrationIndividualsFullyResolvedCollapseAllowNoMigration<-function(popVector,n0multiplierIndividualsList=GenerateN0multiplierIndividuals(popVector,popIntervalsList=GenerateIntervalsFullyResolvedCollapse(popVector,maxK=SetMaxK(popVector)),maxK=SetMaxK(popVector)), maxK=SetMaxK(popVector), verbose=FALSE,file=NULL) {
  migrationArray<-GenerateMigrationIndividualsAllowNoMigration(popVector,n0multiplierIndividualsList, maxK, verbose=verbose, file=NULL)
  migrationIndividualsToKill<-c()
  for (i in sequence(length(migrationArray))) {
     if (RowMax(migrationArray[[i]]$migrationArray[, , 1])==0 && RowMax(migrationArray[[i]]$migrationArray[, , 1])==0) {
        migrationIndividualsToKill<-append( migrationIndividualsToKill,i)
     }
  }
  migrationArray<-migrationArray[-1*migrationIndividualsToKill]
  if (!is.null(file)) {
    save(migrationArray,maxK,popVector,file=file,compress=TRUE)
  }
  return(migrationArray)
}


#this is like generateMigrationIndividuals, but allows some migration rates to be set to 0 with no penalty in terms of number of free parameters
GenerateMigrationIndividualsAllowNoMigration<-function(popVector,n0multiplierIndividualsList=GenerateN0multiplierIndividuals(popVector,popIntervalsList=GenerateIntervals(popVector,maxK),maxK), maxK, verbose=FALSE, file=NULL) {
	migrationIndividualsList<-list()
	for (i in 1:length(n0multiplierIndividualsList)) {
		if (verbose==TRUE) {
			print(paste("doing ",i,"/",length(n0multiplierIndividualsList), " (migration individuals so far = ",length(migrationIndividualsList),")"))
		}
		collapseMatrix<-n0multiplierIndividualsList[[i]]$collapseMatrix
		n0multiplierMap<-n0multiplierIndividualsList[[i]]$n0multiplierMap
		numFinalPops<-dim(collapseMatrix)[1]
		numSteps<-dim(collapseMatrix)[2]
		if ((KCollapseMatrix(collapseMatrix) + KN0multiplierMap(n0multiplierMap) )<=maxK) {
			migrationTemplate<-array(data=NA,dim=c(numFinalPops,numFinalPops,numSteps),dimnames=c("from","to","generation"))
			for (interval in 1:numSteps) {
				for (fromPop in 1:numFinalPops) {
					for (toPop in 1:numFinalPops) {
						if (fromPop!=toPop) {
							if (!is.na(collapseMatrix[fromPop,interval]) && !is.na(collapseMatrix[toPop,interval])) {
								migrationTemplate[fromPop,toPop,interval]<-1
							}
						}
					}
				}
			}
			numPairs=sum(migrationTemplate,na.rm=TRUE)
			mappingMaxKAllowed <- maxK - ( KCollapseMatrix(collapseMatrix) + KN0multiplierMap(n0multiplierMap) )
			evaluateMapping<-TRUE
			if (mappingMaxKAllowed<0) {
				evaluateMapping<-FALSE #since there's no way to do this and not have too many parameters
			}
			thisMapping<-firstcomposition(numPairs)		
			while(evaluateMapping) {
				thisMapping<-c(thisMapping,rep(0,numPairs-length(thisMapping))) #adds trailing zeros
				thisMappingOrig<-thisMapping
				if ((length(which(thisMapping>0)) - 1 ) <=mappingMaxKAllowed) { #note the -1: we take away one free param because we allow zero migration rates
					migrationArray<-migrationTemplate
					whichPositions <- which(migrationArray==1)
					for (positionIndex in 1:length(whichPositions)) {
						position=whichPositions[positionIndex]
						paramPosition<-which(thisMapping>0)[1]
						migrationArray[position]=paramPosition #the position of the first parameter
						thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we have used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
					}
					if (is.na(max(migrationArray,na.rm=TRUE))) {
						migrationIndividualsList[[length(migrationIndividualsList)+1]]<-Migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArray)				
						#print(paste("Just created migration individual ",length(migrationIndividualsList)))
					}
					else {
						for (paramToMakeZero in sequence(max(migrationArray,na.rm=TRUE))) {
							migrationArrayModified<-migrationArray
							migrationArrayModified[which(migrationArrayModified==paramToMakeZero)]<-0
							migrationArrayModified[which(migrationArrayModified>paramToMakeZero)]<-migrationArrayModified[which(migrationArrayModified>paramToMakeZero)]-1
							migrationIndividualsList[[length(migrationIndividualsList)+1]]<-Migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArrayModified)
							#print(paste("Just created migration individual ",length(migrationIndividualsList)))
						}
					}
				}
				if (!islastcomposition(thisMappingOrig,restricted=FALSE)) {
					thisMapping<-nextcomposition(thisMappingOrig,restricted=FALSE)
				}
				else {
					evaluateMapping<-FALSE
				}
			}
		}
	}
  if (!is.null(file)) {
    migrationArray<-migrationIndividualsList
    save(migrationArray,maxK,popVector,file=file,compress=TRUE)
  }
	return(migrationIndividualsList)
}

GenerateExpansionMultiplierIndividuals<-function(popVector,migrationIndividualsList=GenerateMigrationIndividualsAllowNoMigration(popVector,n0multiplierIndividualsList=GenerateN0multiplierIndividuals(popVector,popIntervalsList=GenerateIntervals(popVector,maxK),maxK),maxK),maxK) {
	expansionmultiplierIndividualsList<-list()
	for (i in 1:length(migrationIndividualsList)) {
		expansionmultiplierMapTemplate<-1+0*migrationIndividualsList[[i]]$collapseMatrix  #will have all the populations, all with either NA or 1
		numLineages=sum(expansionmultiplierMapTemplate,na.rm=TRUE)
		possibleMappings<-compositions(numLineages)
		for (mappingIndex in 1:dim(possibleMappings)[2]) {
			thisMapping<-possibleMappings[,mappingIndex]
			if ((length(which(thisMapping>0))+KPopInterval(migrationIndividualsList[[i]]) )<=maxK) { #only do it on those mappings that have not too many free parameters
				expansionmultiplierMap<-expansionmultiplierMapTemplate
				whichPositions <- which(expansionmultiplierMap==1)
				for (positionIndex in 1:length(whichPositions)) {
					position=whichPositions[positionIndex]
					paramPosition<-which(thisMapping>0)[1]
					expansionmultiplierMap[position]=paramPosition #the position of the first parameter
					thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we have used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
				}
				expansionmultiplierIndividualsList[[length(expansionmultiplierIndividualsList)+1]]<-Expansionmultiplierindividual(migrationIndividualsList[[i]]$collapseMatrix, migrationIndividualsList[[i]]$complete, expansionmultiplierMap)
			}
		}
	}
	return(expansionmultiplierIndividualsList)
}


LoadMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation="/usr/local/bin/ms") {
	msCallInfo<-CreateMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
	geneTrees<-system(paste(msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';'"),intern=TRUE)
	return(geneTrees)
}

SaveMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation="/usr/local/bin/ms",file="sim.tre") {
	msCallInfo<-CreateMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
	returnCode<-system(paste(msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';' >",file),intern=FALSE)
	return(returnCode)
}

PipeMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, debug=FALSE,print.ms.string=FALSE,ncores=1) {
	msCallInfo<-CreateMSstringSpecific(popVector,migrationIndividual,parameterVector,ceiling(nTrees/ncores))
	
  if(print.ms.string) {
    print(msCallInfo) 
  }
  if(debug) {
    print("parameterVector")
    print(parameterVector)
  }
	unresolvedFlag<-"-u"
	if (unresolvedTest==FALSE) {
		unresolvedFlag<-""
	}
	outputstring<-paste(msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';' | perl ",compareLocation, unresolvedFlag, paste("-a",assign,sep=""), paste("-o",observed,sep=""), sep=" ")
	if (debug) {
		print(outputstring) 
	}
	#if (ncores==1) {
		outputVector<-system(outputstring,intern=TRUE)
		return(outputVector)
	#} else {
	#	wrapOutput<-function(x,outputstring) {
	#		as.numeric(system(outputstring,intern=TRUE))
	#	}
	#	outputVector<-apply(simplify2array(mclapply(sequence(ncores),wrapOutput,outputstring=outputstring,mc.cores=ncores)),1,sum)
	#	return(outputVector)
	#}
}


TaxaToRetain<-function(assignFrame,nIndividualsDesired,minPerPop=1,attemptsCutoff=100000,finalPopVector=NULL) {
	#if "finalPopVector" has not been specified, then iteratively randomly sample "nIndividualsDesired" from the 
	#entire dataset "attemptsCutoff" times until at least "minPerPop" individuals are represented per population 
	#within "toRetain"	
	if(is.null(finalPopVector)) {
		samplesGood<-FALSE
		toRetain<-c()
		attempts<-0
		while (samplesGood!=TRUE) {
			if (attempts>attemptsCutoff) {
				stop(paste("No random sample met the criteria after",attempts,"attempts"))
			}		
			attempts<-attempts+1
			toRetain<-sample(dim(assignFrame)[1],nIndividualsDesired,replace=FALSE)
			retainedPops<-(assignFrame[,1])[toRetain]
			samplesGood<-FALSE
			if (nlevels(retainedPops)==nlevels(assignFrame[,1])) {
				samplesGood<-TRUE
				for (i in 1:nlevels(assignFrame[,1])){
					if (length(which(retainedPops==(levels(assignFrame[,1]))[i]))<minPerPop) {
						samplesGood<-FALSE
					}
				}
			}
		}
	} else {
		#if finalPopVector has been specified, then randomly sample according to the population-specific sample sizes
		#indicated within finalPopVector
		nIndividualsDesired<-sum(finalPopVector)
		toRetain <- array()
		#for each population...
		for(numLevel in 1:length(levels(assignFrame[,1])))
		{
			thisPopSamples <- assignFrame[which(assignFrame[,1] == unique(assignFrame[,1])[numLevel]),][[2]]
			#...if the there is more than one individual in the population...
			if (length(thisPopSamples) > 1)
			{
				#...randomly sample "finalPopVector[numLevel]" individuals and add them to "toRetain"
				toRetainTemp <- sample(thisPopSamples,finalPopVector[numLevel],replace=FALSE)
				toRetain <- c(toRetain,toRetainTemp)[!is.na(c(toRetain,toRetainTemp))]
			} else {
				#But if there is only one individual in the population (and one sample is to be drawn)...
				if(finalPopVector[numLevel] == 1)
				{
					#Then add that sample to "toRetain"
					toRetain <- c(toRetain,thisPopSamples)[!is.na(c(toRetain,thisPopSamples))]
				}
			}
		}
	}
return(toRetain[sort.list(toRetain)])
}

TaxaToDrop<-function(assignFrame,taxaRetained) {
	allTaxa<-c(1:dim(assignFrame)[1])
	taxaToDrop<-allTaxa[-taxaRetained]
	return(taxaToDrop)
}

function(assignFrame,phy,nIndividualsDesired,nSamplesDesired,minPerPop=1,finalPopVector=NULL,subsamplingPath=NULL,
	observedSubsampleFile="observed",assignSubsampleFile="assign",outgroupPrune=TRUE) {
	retainedTaxaMatrix<-matrix(NA,nrow=nSamplesDesired,ncol=nIndividualsDesired)
	#Iteratively subsamples using TaxaToRetain and then prunes the tree and assingment file accordingly. If outgroupPrune is
	#TRUE, then the last individual in the tree is removed.
	if(class(phy)!="multiPhylo") {
		phy<-c(phy)
	}
	for (rep in 1:nSamplesDesired) {
		keepTaxa<-TaxaToRetain(assignFrame,nIndividualsDesired,minPerPop,attemptsCutoff=100000,finalPopVector)
		retainedTaxaMatrix[rep,]<-keepTaxa
		prunedAF<-PrunedAssignFrame(assignFrame,keepTaxa)
		physamp<-phy
		delTaxa<-TaxaToDrop(assignFrame,keepTaxa)
		for (tree in 1:length(phy)) {
			newphy<-drop.tip(phy[[tree]],as.character(delTaxa))
			for (tipIndex in 1:length(newphy$tip.label)) {
				old.label<-newphy$tip.label[tipIndex]
				new.label<-as.character(which(prunedAF[,2]==old.label))
				newphy$tip.label[tipIndex]<-new.label
			}
			if(outgroupPrune==TRUE){
				physamp[[1]]<-drop.tip(newphy,length(keepTaxa))
			}else{			
			physamp[[tree]]<-newphy
			}
		}
		prunedAF[,2]<-as.character(c(1:length(prunedAF[,2])))
		if(outgroupPrune==TRUE){
			prunedAF<-prunedAF[-length(prunedAF[,2]),]
		}
		if(is.null(subsamplingPath)) {
			write.tree(physamp,file=paste(observedSubsampleFile,rep,".tre",sep=""))
			write.table(prunedAF,file=paste(assignSubsampleFile,rep,".txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
		}
		else {
			write.tree(physamp,file=paste(subsamplingPath,observedSubsampleFile,rep,".tre",sep=""))
			write.table(prunedAF,file=paste(subsamplingPath,assignSubsampleFile,rep,".txt",sep=""),
				quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
			
		}
	}
	return(retainedTaxaMatrix)
}


PrepAssignFrame <- function(phy,assignFrame) {
#Function for prepping the assignFrame for subsampling
#This function inputs 1) an assignment file that includes all samples pooled from across loci and 2) trees for each locus and outputs
#assignment files tailored to each locus and trees that are relabeled with numbers to match the assignment files. Note that only tree tips
#that are included in the inputted assignment file will be retained in the outputs.
 	
	tips.vec <- phy$tip.label #make an array of the included individuals
	phy$tip.label <- as.character(c(1:length(phy$tip.label))) #change tip labels to consecutive numbers			
	assign.vec <- as.character(assignFrame[,1]) #make an array of all individuals in the dataset
	match.vec <- data.frame()
	for(match in 1:length(tips.vec)) #for each tip in the tree, in order (so that tip labels can be changed to numbers and keep the same order)
	{
		match.temp = NULL
		for(match1 in 1:length(assign.vec)) #go through each individual in the assignment spreadsheet and ask...
		{
			a.match <- tips.vec[match] %in% assign.vec[match1] #is the current tip individual the same as the current individual in the spreadsheet?
			if(a.match=="TRUE") #if so,
			{
				match.temp <- assignFrame[match1,] #remember the assignment
			}
		}
		if(!is.null(match.temp)){
			match.vec <- rbind(match.vec,match.temp) #add this assignment to a locus-specific assignment spreadsheet
		} else {
			phy<-drop.tip(phy,as.character(match)) #unless there was no matching sample in the assingment file, in which case, drop this tip from the tree
			cat("Warning: Some trees contain tip names not included in the inputted assignment file. These tips will not be subsampled.\n",
				file=paste(basePath,"WARNINGS.txt",sep=""),append=TRUE)
		}
	}
	assignFrame<-rbind(data.frame(match.vec$popLabel,c(1:length(phy$tip.label)))) #convert this locus-specific spreadsheet into assignFrame
	colnames(assignFrame) <- c("popLabel","indivTotal")	
	assignFrame <- assignFrame[order(assignFrame$popLabel),] #order new assignFrame by population
	
	#Now re-name the tree tips so they match the assignment file order (tips belonging to population A are numbered first, B second, C third, etc)
	for(changetips in 1:length(phy$tip.label))
	{
		phy$tip.label[assignFrame$indivTotal[changetips]] <- as.numeric(changetips)
	}

	#Now, make new assignFrame with the new tip labels listed in order
	assignFrame <- data.frame(assignFrame[,1],c(1:length(phy$tip.label)))
	colnames(assignFrame) <- c("popLabel","indivTotal")	
	return(list(assignFrame,phy)) #return both the assignFrame and re-labeled trees
}	

PrunedPopVector<-function(assignFrame,taxaRetained) {
	assignFrameLabels<-assignFrame[taxaRetained,1]
	popVector<-rep(NA,nlevels(assignFrameLabels))
	for (i in 1:nlevels(assignFrameLabels)) {
		popVector[i]<-length(which(assignFrameLabels==levels(assignFrameLabels)[i])) 
	}
	return(popVector)
}

PrunedAssignFrame<-function(assignFrame,taxaRetained) {
	return(assignFrame[taxaRetained,])
}

#idea is that you scale the prob of missing by the number of possible trees
ConvertOutputVectorToLikelihood<-function(outputVector,nTrees,probOfMissing=(1/howmanytrees(sum(popVector)))) {
	outputVector<-as.numeric(outputVector)
	outputVector[which(outputVector==0)]<-probOfMissing
	outputVector<-outputVector/nTrees
	outputVector<-log(outputVector)
	lnL<-sum(outputVector)
	return(lnL)
}

#Uses ln lik
CombineSubsampleLikelihoods<-function(likelihoodVector,nIndividualsDesired,orig.popVector) {
	originalSize<-sum(orig.popVector)
	originalNumberTrees<-howmanytrees(originalSize)
	prunedNumberTrees<-howmanytrees(nIndividualsDesired)
	sumLnL<-sum(likelihoodVector)
	meanLnL<-mean(likelihoodVector)
	medianLnL<-median(likelihoodVector)
#prob(big tree)=prob(small tree)*(#small trees / #big trees)
#ln(prob(big tree) = ln( prob(small tree)*(#small trees / #big trees) )
#ln(prob(big tree) = ln( prob(small tree) ) + ln (#small trees)  - ln (#big trees)
	rescaledLikelihoodVector<-likelihoodVector + log(prunedNumberTrees) - log(originalNumberTrees)
	rescaledSumLnL<-sum(rescaledLikelihoodVector)
	rescaledMeanLnL<-mean(rescaledLikelihoodVector)
	rescaledMedianLnL<-median(rescaledLikelihoodVector)
	return(list(sumLnL,meanLnL,medianLnL,rescaledSumLnL,rescaledMeanLnL,rescaledMedianLnL))
}

GenerateMigrationArrayMap<-function(migrationArray) {
  migrationArrayMap<-data.frame(matrix(c(1,1,1,1,1,1,1),nrow=1)) #uses the first entry in migrationArray
  print(migrationArrayMap)
  for (model in 2:length(migrationArray)) {
    newRow<-c(model,rep(NA,6))
    for (comparison in sequence(dim(migrationArrayMap)[1])) {
      comparisonRow<-as.integer(migrationArrayMap[comparison,])
     if(identical(migrationArray[[comparisonRow[3] ]]$collapseMatrix, migrationArray[[model]]$collapseMatrix)) {
       newRow[2:3]<-comparisonRow[2:3]
     }
     if(identical(migrationArray[[comparisonRow[5] ]]$n0multiplierMap, migrationArray[[model]]$n0multiplierMap)) {
       newRow[4:5]<-comparisonRow[4:5]
     }
     if(identical(migrationArray[[comparisonRow[7] ]]$migrationArray, migrationArray[[model]]$migrationArray)) {
       newRow[6:7]<-comparisonRow[6:7]
     }
    }
    if(is.na(newRow[2])) {
      newRow[2:3]<-c(1+max(migrationArrayMap[,2]),model)
    }
    if(is.na(newRow[4])) {
      newRow[4:5]<-c(1+max(migrationArrayMap[,4]),model)
    }
   if(is.na(newRow[6])) {
      newRow[6:7]<-c(1+max(migrationArrayMap[,6]),model)
    }
    migrationArrayMap<-rbind(migrationArrayMap,newRow)
    print(newRow)
  }
  names(migrationArrayMap)<-c("model","collapseMatrix.number","collapseMatrix.parent","n0multiplierMap.number","n0multiplierMap.parent","migrationArray.number","migrationArray.parent")
  return(migrationArrayMap)
}

#only works on island models (no collapse)
ReturnSymmetricMigrationOnly<-function(migrationArray) {
  newMigrationArray<-list()
  storedObjectString<-c()
  for (i in sequence(length(migrationArray))) {
    migrationIndividual<-migrationArray[[i]]
    if(max(migrationIndividual$collapseMatrix,na.rm=TRUE)==0) {
      migration.main<-migrationIndividual$migrationArray[ , , 1] 
      migration.upper<-migration.main[upper.tri(migration.main)]
      if(length(c(min(migration.upper):max(migration.upper)))==length(unique(migration.upper))) { #that is, we don't skip any parameters: no 0, 1, 3
        for (row.index in 2:dim(migration.main)[1]) {
          for (col.index in 1:(dim(migration.main)[2]-1)) {
            if (row.index>col.index) {
              migration.main[row.index,col.index]=migration.main[col.index,row.index] 
            }
          }
        }
        migrationIndividual$migrationArray[, , 1]<-migration.main
        focalString<-paste(as.character(dput(migrationIndividual)),collapse="_")
        if(MatrixContainsAllValues(migration.main,1)) {
		   if(sum(grepl(focalString,storedObjectString))==0) {
			 if (min(ColMax(migration.main))>0) { #false means that the maximum value for one column is zero, indicating no inflow into that pop
				if (min(RowMax(migration.main))>0) { #false means that the maximum value for the row is zero, indicating no outflow from that pop (and note that we're supposed to be symmetric, so the colMax should have checked this, but redundancy for safety is good).
					storedObjectString<-c(storedObjectString,focalString)
					newMigrationArray[[1+length(newMigrationArray)]]<-migrationIndividual
				}
          	}
          }
        }
      }
    }
  }
  return(newMigrationArray)
}

RunSeqConverter <- function(seqConvPath,inFile,outputFormat,inputFormat){
	#Function for running seqConverter to convert nexus files to phylip (requires that the seqConverter perl script
	#be in your designated path)
	seqConvPath=seqConvPath
	inFile=inFile
	outputFormat=outputFormat
	inputFormat=inputFormat
	
	systemCall1 <-system(paste(seqConvPath,"seqConverter.pl -d",inFile," -o",outputFormat," -i",inputFormat,sep=""))
	return(systemCall1)
}

RunRaxml <- function(raxmlPath,raxmlVersion,inputPath,inputFile,mutationModel,outgroup,iterations,seed){
	#Function for producing input string for RAxML and calling up the program (requires that RAxML be in your
	#designated path)
	raxmlPath=raxmlPath
	raxmlVersion=raxmlVersion
	inputPath=inputPath
	inputFile=inputFile
	mutationModel=mutationModel
	outgroup=outgroup
	iterations=iterations
	seed=seed
	outputFile <- paste(inputFile,".tre",sep="")
		
	systemCall2 <- system(paste(raxmlPath,raxmlVersion," -w ",inputPath," -s ",inputPath,inputFile," -n ",outputFile," -m ",mutationModel,
		" -o ",outgroup," -f d -N ",iterations," -p ",seed,sep=""))
	return(systemCall2)
}
