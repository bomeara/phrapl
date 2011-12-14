library(partitions)
library(ape)
library(rgenoud)
library(optimx)
library(parallel)


colMax<-function(x,na.rm=TRUE) {
	maxVal=rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		maxVal[i]<-max(x[,i],na.rm=na.rm)
	}
	return(maxVal)
}

colMin<-function(x,na.rm=TRUE) {
	minVal=rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		minVal[i]<-min(x[,i],na.rm=na.rm)
	}
	return(minVal)
}

rowMax<-function(x,na.rm=TRUE) {
  maxVal=rep(NA,dim(x)[1])
	for (i in 1:dim(x)[1]) {
		maxVal[i]<-max(x[i,],na.rm=na.rm)
	}
	return(maxVal)
}

rowMin<-function(x,na.rm=TRUE) {
	minVal=rep(NA,dim(x)[1])
	for (i in 1:dim(x)[1]) {
		minVal[i]<-min(x[i,],na.rm=na.rm)
	}
	return(minVal)
}


colCountIf<-function(x,val) {
	countVec<-rep(NA,dim(x)[2])
	for (i in 1:dim(x)[2]) {
		countVec[i]=length(which(x[,i]==val))
	}
	return(countVec)
}

colCountLast<-function(x,val) {
	return(length(which(x[,dim(x)[2]]==val)))
}

setMaxK<-function(popVector,nTrees=1,taxaPerParam=5,force=FALSE) {
  if (exists("maxK") && !force) {
    print(paste("note: setMaxK called, but maxK already exists, so using the existing one of ",maxK,". To change this behavior, give setMaxK the option force=TRUE",sep=""))
    return(maxK) 
  }
  else {
  	maxK<-max(1,floor(sum(popVector)/(nTrees*taxaPerParam))) 
  	return(maxK)
  }
}

popinterval<-function(collapseMatrix,complete=FALSE) {
#collapse matrix has populations as rows and each generation going back in time as columns
	pInt<-list(collapseMatrix=collapseMatrix,complete=complete)
	class(pInt)="popinterval"
	if (max(collapseMatrix,na.rm=TRUE)==0) {
		pInt$complete=TRUE #if the model is one of no collapse, nothing further to do
	}
	return(pInt)
}

n0multiplierindividual<-function(collapseMatrix,complete=FALSE,n0multiplierMap) {
	tI<-list(collapseMatrix=collapseMatrix,complete=complete,n0multiplierMap=n0multiplierMap)
	class(tI)<-"n0multiplierindividual"
	return(tI)
}

migrationindividual<-function(collapseMatrix,complete=FALSE,n0multiplierMap,migrationArray) {
	tI<-list(collapseMatrix=collapseMatrix,complete=complete,n0multiplierMap=n0multiplierMap,migrationArray=migrationArray)
	class(tI)<-"migrationindividual"
	return(tI)
}

msIndividualParameters<-function(migrationIndividual) {
	collapseMatrix<-migrationIndividual$collapseMatrix
	complete<-migrationIndividual$complete
	n0multiplierMap<-migrationIndividual$n0multiplierMap
	migrationArray<-migrationIndividual$migrationArray
	parameterList<-c()
	if (max(collapseMatrix,na.rm=TRUE)>0) {
		for (i in sequence(kCollapseMatrix(collapseMatrix))) {
			parameterList<-append(parameterList,paste("collapse_",i,sep="")) 
		}
	}
	for (i in sequence(max(n0multiplierMap,na.rm=TRUE))) {
		parameterList<-append(parameterList,paste("n0multiplier_",i,sep="")) 
	}
	for (i in sequence(max(migrationArray,na.rm=TRUE))) {
		parameterList<-append(parameterList,paste("migration_",i,sep="")) 
	}
	return(parameterList)
}

kAll<-function(migrationIndividual) {
  return(length(msIndividualParameters(migrationIndividual))) 
}

createMSstringSpecific<-function(popVector,migrationIndividual,parameterVector,nTrees=1) {
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

createMSstringGeneric<-function(popVector,migrationIndividual,parameterVector,nTrees=1) {
	return(createMSstringSpecific(popVector,migrationIndividual,msIndividualParameters(migrationIndividual),nTrees)) 
}

returnIncompletes<-function(popIntervalsList) {
	incompleteElements<-c()
	for (i in 1:length(popIntervalsList)) {
		if (!popIntervalsList[[i]]$complete) {
			incompleteElements<-append(incompleteElements,i)
		}
	}
	return(incompleteElements)
}

updateCompletes<-function(popIntervalsList) {
	for (i in 1:length(popIntervalsList)) {
		if (colCountLast(popIntervalsList[[i]]$collapseMatrix,1)==0) {
			popIntervalsList[[i]]$complete<-TRUE #end up with having no more collapses
		}
		else if (colCountLast(popIntervalsList[[i]]$collapseMatrix,0)==0) {
			popIntervalsList[[i]]$complete<-TRUE #end up as one population
		}
	}
	return(popIntervalsList)	
}

completeIntervals<-function(popIntervalsList) {
	while(length(returnIncompletes(popIntervalsList))>0) {
		incompleteElements=returnIncompletes(popIntervalsList)
		newPopIntervalsList<-popIntervalsList[-1*incompleteElements] #has only complete elements
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
				rawIntervals<-rawIntervals[,colMax(rawIntervals)<=1] #we are okay with having all zeros: no population collapse
				rawIntervals<-rawIntervals[,colCountIf(rawIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
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
				newPopIntervalsList[[length(newPopIntervalsList)+1]]<-popinterval(cbind(currentParent$collapseMatrix,matrix(rawIntervalsRescaled[,k],ncol=1)))
			}
		}
		popIntervalsList<-newPopIntervalsList
		popIntervalsList<-updateCompletes(popIntervalsList)
	}
	return(popIntervalsList)
}


generateIntervals<-function(popVector,maxK) {
#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
#maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
	nPop <- length(popVector)
	firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	firstIntervals<-firstIntervals[,colMax(firstIntervals)<=1] #we are okay with having all zeros: no population collapse
	firstIntervals<-firstIntervals[,colCountIf(firstIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
	popIntervalsList<-list()
	for (i in 1:dim(firstIntervals)[2]) {
		popIntervalsList[[i]]<-popinterval(as.matrix(firstIntervals[,i],ncol=1))
	}
	popIntervalsList<-completeIntervals(updateCompletes(popIntervalsList))
	return(popIntervalsList)
}

generateIntervalsNoCollapse<-function(popVector,maxK) {
#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
#maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
  nPop <- length(popVector)
	firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	firstIntervals<-firstIntervals[,colMax(firstIntervals)==0] #we only want ones with no collapse (like an island model)
	popIntervalsList<-list()
	popIntervalsList[[1]]<-popinterval(as.matrix(firstIntervals,ncol=1))
	popIntervalsList<-completeIntervals(updateCompletes(popIntervalsList))
	return(popIntervalsList)
}

generateIntervalsFullyResolvedCollapse<-function(popVector,maxK) {
#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
#maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
  nPop <- length(popVector)
	firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	firstIntervals<-firstIntervals[,colMax(firstIntervals)<=1] #we are okay with having all zeros: no population collapse
	firstIntervals<-firstIntervals[,colCountIf(firstIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
	popIntervalsList<-list()
	for (i in 1:dim(firstIntervals)[2]) {
		popIntervalsList[[i]]<-popinterval(as.matrix(firstIntervals[,i],ncol=1))
	}
	popIntervalsList<-completeIntervals(updateCompletes(popIntervalsList))
  intervalsToDelete<-c()
  for (i in sequence(length(popIntervalsList))) {
    focalCollapseMatrix<-popIntervalsList[[i]]$collapseMatrix
    if ((min(colMax(focalCollapseMatrix))==0) || (dim(focalCollapseMatrix)[1]!=(1+dim(focalCollapseMatrix)[2]))) {
       intervalsToDelete<-append(intervalsToDelete,i)
    }
  }
  popIntervalsList<-popIntervalsList[-intervalsToDelete]
	return(popIntervalsList)
}

kPopInterval<-function(popInterval) {
#returns the number of free parameters needed for that interval object. For example, having everything collapse in one step requires one param (the tMRCA). Having one collapse then a second requires two. Having no collapse requires 0
	maxVector<-colMax(popInterval$collapseMatrix)
	return(length(which(maxVector>0)))
}

kCollapseMatrix<-function(collapseMatrix) {
#returns the number of free parameters needed for that interval object. For example, having everything collapse in one step requires one param (the tMRCA). Having one collapse then a second requires two. Having no collapse requires 0
	maxVector<-colMax(collapseMatrix)
	return(length(which(maxVector>0)))
}

#the basic idea here is that at each population in each time interval there is a n0multiplier. These can all be set to the same value, allowed to vary, or assigned in clumps (i.e., pops 1, 3, and 6 have the same n0multiplier value)
#this generates all such mappings, subject to staying within the maximum number of free parameters
generateN0multiplierIndividuals<-function(popVector,popIntervalsList=generateIntervals(popVector,maxK=setMaxK(popVector)),maxK=setMaxK(popVector)) {
	n0multiplierIndividualsList<-list()
	for (i in 1:length(popIntervalsList)) {
		n0multiplierMapTemplate<-1+0*popIntervalsList[[i]]$collapseMatrix  #will have all the populations, all with either NA or 1
		numLineages=sum(n0multiplierMapTemplate,na.rm=TRUE)
		possibleMappings<-compositions(numLineages)
		for (mappingIndex in 1:dim(possibleMappings)[2]) {
			thisMapping<-possibleMappings[,mappingIndex]
			if ((length(which(thisMapping>0))+kPopInterval(popIntervalsList[[i]]) )<=maxK) { #only do it on those mappings that have not too many free parameters
				n0multiplierMap<-n0multiplierMapTemplate
				whichPositions <- which(n0multiplierMap==1)
				for (positionIndex in 1:length(whichPositions)) {
					position=whichPositions[positionIndex]
					paramPosition<-which(thisMapping>0)[1]
					n0multiplierMap[position]=paramPosition #the position of the first parameter
					thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we have used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
				}
				n0multiplierIndividualsList[[length(n0multiplierIndividualsList)+1]]<-n0multiplierindividual(popIntervalsList[[i]]$collapseMatrix, popIntervalsList[[i]]$complete, n0multiplierMap)
			}
		}
	}
	return(n0multiplierIndividualsList)
}

kN0multiplierMap<-function(kN0multiplierMap) {
	return(max(n0multiplierMap,na.rm=TRUE)) 
}

#now we will generate all possible assignments of pairwise migration. Again, we want to keep the total number of free parameters (times, n0multipliers, migration rates) under our chosen max
#allow a model where migrations change anywhere along branch, or only at coalescent nodes? The problem with the latter is that you cannott fit some reasonable models: i.e., two populations persisting through time. Problem with the former is parameter space
generateMigrationIndividuals<-function(popVector,n0multiplierIndividualsList=generateN0multiplierIndividuals(popVector,popIntervalsList=generateIntervals(popVector,maxK=setMaxK(popVector)),maxK=setMaxK(popVector)), maxK=setMaxK(popVector), verbose=FALSE) {
	migrationIndividualsList<-list()
	for (i in sequence(n0multiplierIndividualsList)) {
		if (verbose==TRUE) {
			print(paste("doing ",i,"/",length(n0multiplierIndividualsList)))
		}
		collapseMatrix<-n0multiplierIndividualsList[[i]]$collapseMatrix
		n0multiplierMap<-n0multiplierIndividualsList[[i]]$n0multiplierMap
		numFinalPops<-dim(collapseMatrix)[1]
		numSteps<-dim(collapseMatrix)[2]
		if ((kCollapseMatrix(collapseMatrix) + kN0multiplierMap(n0multiplierMap) )<=maxK) {
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
			mappingMaxKAllowed <- maxK - ( kCollapseMatrix(collapseMatrix) + kN0multiplierMap(n0multiplierMap) )
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
					migrationIndividualsList[[length(migrationIndividualsList)+1]]<-migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArray)
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
generateMigrationIndividualsNoCollapseAllowNoMigration<-function(popVector,n0multiplierIndividualsList=generateN0multiplierIndividuals(popVector,popIntervalsList=generateIntervalsNoCollapse(popVector,maxK=setMaxK(popVector)),maxK=setMaxK(popVector)), maxK=setMaxK(popVector), verbose=FALSE,file=NULL) {
  migrationArray<-generateMigrationIndividualsAllowNoMigration(popVector,n0multiplierIndividualsList, maxK, verbose=verbose, file=NULL)
  migrationIndividualsToKill<-c()
  for (i in sequence(length(migrationArray))) {
     if (rowMax(migrationArray[[i]]$migrationArray[, , 1])==0 && rowMax(migrationArray[[i]]$migrationArray[, , 1])==0) {
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
generateMigrationIndividualsFullyResolvedCollapseAllowNoMigration<-function(popVector,n0multiplierIndividualsList=generateN0multiplierIndividuals(popVector,popIntervalsList=generateIntervalsFullyResolvedCollapse(popVector,maxK=setMaxK(popVector)),maxK=setMaxK(popVector)), maxK=setMaxK(popVector), verbose=FALSE,file=NULL) {
  migrationArray<-generateMigrationIndividualsAllowNoMigration(popVector,n0multiplierIndividualsList, maxK, verbose=verbose, file=NULL)
  migrationIndividualsToKill<-c()
  for (i in sequence(length(migrationArray))) {
     if (rowMax(migrationArray[[i]]$migrationArray[, , 1])==0 && rowMax(migrationArray[[i]]$migrationArray[, , 1])==0) {
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
generateMigrationIndividualsAllowNoMigration<-function(popVector,n0multiplierIndividualsList=generateN0multiplierIndividuals(popVector,popIntervalsList=generateIntervals(popVector,maxK),maxK), maxK, verbose=FALSE, file=NULL) {
	migrationIndividualsList<-list()
	for (i in 1:length(n0multiplierIndividualsList)) {
		if (verbose==TRUE) {
			print(paste("doing ",i,"/",length(n0multiplierIndividualsList), " (migration individuals so far = ",length(migrationIndividualsList),")"))
		}
		collapseMatrix<-n0multiplierIndividualsList[[i]]$collapseMatrix
		n0multiplierMap<-n0multiplierIndividualsList[[i]]$n0multiplierMap
		numFinalPops<-dim(collapseMatrix)[1]
		numSteps<-dim(collapseMatrix)[2]
		if ((kCollapseMatrix(collapseMatrix) + kN0multiplierMap(n0multiplierMap) )<=maxK) {
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
			mappingMaxKAllowed <- maxK - ( kCollapseMatrix(collapseMatrix) + kN0multiplierMap(n0multiplierMap) )
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
						migrationIndividualsList[[length(migrationIndividualsList)+1]]<-migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArray)				
						#print(paste("Just created migration individual ",length(migrationIndividualsList)))
					}
					else {
						for (paramToMakeZero in sequence(max(migrationArray,na.rm=TRUE))) {
							migrationArrayModified<-migrationArray
							migrationArrayModified[which(migrationArrayModified==paramToMakeZero)]<-0
							migrationArrayModified[which(migrationArrayModified>paramToMakeZero)]<-migrationArrayModified[which(migrationArrayModified>paramToMakeZero)]-1
							migrationIndividualsList[[length(migrationIndividualsList)+1]]<-migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArrayModified)
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

generateExpansionMultiplierIndividuals<-function(popVector,migrationIndividualsList=generateMigrationIndividualsAllowNoMigration(popVector,n0multiplierIndividualsList=generateN0multiplierIndividuals(popVector,popIntervalsList=generateIntervals(popVector,maxK),maxK),maxK),maxK) {
	expansionmultiplierIndividualsList<-list()
	for (i in 1:length(migrationIndividualsList)) {
		expansionmultiplierMapTemplate<-1+0*migrationIndividualsList[[i]]$collapseMatrix  #will have all the populations, all with either NA or 1
		numLineages=sum(expansionmultiplierMapTemplate,na.rm=TRUE)
		possibleMappings<-compositions(numLineages)
		for (mappingIndex in 1:dim(possibleMappings)[2]) {
			thisMapping<-possibleMappings[,mappingIndex]
			if ((length(which(thisMapping>0))+kPopInterval(migrationIndividualsList[[i]]) )<=maxK) { #only do it on those mappings that have not too many free parameters
				expansionmultiplierMap<-expansionmultiplierMapTemplate
				whichPositions <- which(expansionmultiplierMap==1)
				for (positionIndex in 1:length(whichPositions)) {
					position=whichPositions[positionIndex]
					paramPosition<-which(thisMapping>0)[1]
					expansionmultiplierMap[position]=paramPosition #the position of the first parameter
					thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we have used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
				}
				expansionmultiplierIndividualsList[[length(expansionmultiplierIndividualsList)+1]]<-expansionmultiplierindividual(migrationIndividualsList[[i]]$collapseMatrix, migrationIndividualsList[[i]]$complete, expansionmultiplierMap)
			}
		}
	}
	return(expansionmultiplierIndividualsList)
}


loadMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation="/usr/local/bin/ms") {
	msCallInfo<-createMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
	geneTrees<-system(paste(msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';'"),intern=TRUE)
	return(geneTrees)
}

saveMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation="/usr/local/bin/ms",file="sim.tre") {
	msCallInfo<-createMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
	returnCode<-system(paste(msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';' >",file),intern=FALSE)
	return(returnCode)
}

pipeMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, debug=FALSE,print.ms.string=FALSE) {
	msCallInfo<-createMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
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
	outputVector<-system(outputstring,intern=TRUE)
	return(outputVector)
}



createAssignment<-function(popVector,file="assign.txt") {
	alphabet<-strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ","")[[1]]
	assignFrame<-data.frame()
	indivTotal<-0
	for(popIndex in 1:length(popVector)) {
		popLabel<-alphabet[popIndex]
		for(indivIndex in 1:popVector[popIndex]) {
			indivTotal<-indivTotal+1
			if(indivTotal==1) {
				assignFrame<-data.frame(popLabel,indivTotal) 
			}
			else {
				assignFrame<-rbind(assignFrame,data.frame(popLabel, indivTotal))
			}
		}
	}
	write.table(assignFrame,file=file,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}

createAssignment.df<-function(popVector) {
	alphabet<-strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ","")[[1]]
	assignFrame<-data.frame()
	indivTotal<-0
	for(popIndex in 1:length(popVector)) {
		popLabel<-alphabet[popIndex]
		for(indivIndex in 1:popVector[popIndex]) {
			indivTotal<-indivTotal+1
			if(indivTotal==1) {
				assignFrame<-data.frame(popLabel,indivTotal) 
			}
			else {
				assignFrame<-rbind(assignFrame,data.frame(popLabel, indivTotal))
			}
		}
	}
	return(assignFrame)
}

taxaToRetain<-function(assignFrame,nIndividualsDesired,minPerPop=1,attemptsCutoff=100000) {
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
	return(toRetain[sort.list(toRetain)])
}

taxaToDrop<-function(assignFrame,taxaRetained) {
	allTaxa<-c(1:dim(assignFrame)[1])
	taxaToDrop<-allTaxa[-taxaRetained]
	return(taxaToDrop)
}

prepSubsampling<-function(assignFrame,phy, nIndividualsDesired,nSamplesDesired,minPerPop=1) {
	retainedTaxaMatrix<-matrix(NA,nrow=nSamplesDesired,ncol=nIndividualsDesired)
	for (rep in 1:nSamplesDesired) {
		keepTaxa<-taxaToRetain(assignFrame,nIndividualsDesired,minPerPop)
		retainedTaxaMatrix[rep,]<-keepTaxa
		prunedAF<-prunedAssignFrame(assignFrame,keepTaxa)
		physamp<-phy
		delTaxa<-taxaToDrop(assignFrame,keepTaxa)
		for (tree in 1:length(phy)) {
			newphy<-drop.tip(phy[[tree]],as.character(delTaxa))
			for (tipIndex in 1:length(newphy$tip.label)) {
				old.label<-newphy$tip.label[tipIndex]
# print(paste("old.label = ",old.label))
# print(prunedAF[,2])
				new.label<-as.character(which(prunedAF[,2]==old.label))
				newphy$tip.label[tipIndex]<-new.label
			}
			physamp[[tree]]<-newphy
		}
		prunedAF[,2]<-as.character(c(1:length(prunedAF[,2])))
		write.tree(physamp,file=paste("obs",rep,".tre",sep=""))
		write.table(prunedAF,file=paste("assign",rep,".txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	}
	return(retainedTaxaMatrix)
}

prunedPopVector<-function(assignFrame,taxaRetained) {
	assignFrameLabels<-assignFrame[taxaRetained,1]
	popVector<-rep(NA,nlevels(assignFrameLabels))
	for (i in 1:nlevels(assignFrameLabels)) {
		popVector[i]<-length(which(assignFrameLabels==levels(assignFrameLabels)[i])) 
	}
	return(popVector)
}

prunedAssignFrame<-function(assignFrame,taxaRetained) {
	return(assignFrame[taxaRetained,])
}

convertOutputVectorToLikelihood<-function(outputVector,nTrees) {
	outputVector<-as.numeric(outputVector)
	outputVector[which(outputVector==0)]<-0.1
	outputVector<-outputVector/nTrees
	outputVector<-log(outputVector)
	lnL<-sum(outputVector)
	return(lnL)
}

#Uses ln lik
combineSubsampleLikelihoods<-function(likelihoodVector,nIndividualsDesired,orig.popVector) {
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

generateMigrationArrayMap<-function(migrationArray) {
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

returnModel<-function(p,migrationArrayMap) {
   prunedResults<-subset(migrationArrayMap, migrationArrayMap$collapseMatrix.number==p[1])
   prunedResults<-subset(prunedResults, prunedResults$n0multiplierMap.number==p[2])
   prunedResults<-subset(prunedResults, prunedResults$migrationArray.number==p[3])
   if(dim(prunedResults)[1]==1) {
      return(prunedResults$model) 
   }
   else  {
      return(NA) 
   }
}

#maxParameterValue prevents MS from going nuts (not finishing) with really high migration or other rates
returnAIC<-function(par,popVector,migrationIndividual,nTrees=1,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, print.results=FALSE, print.ms.string=FALSE, debug=FALSE, badAIC=100000000000000, maxParameterValue=1000,...) {
  parameterVector<-exp(par)
  names(parameterVector)<-msIndividualParameters(migrationIndividual)
  if(max(parameterVector)>maxParameterValue) {
  	return(badAIC)
  }
  if(debug) {
    print(parameterVector) 
  }
  lnLValue<-convertOutputVectorToLikelihood(pipeMS(popVector=popVector,migrationIndividual=migrationIndividual,parameterVector=parameterVector,nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest,print.ms.string=print.ms.string, debug=debug),nTrees=nTrees)
  AICValue<-2*(-lnLValue + kAll(migrationIndividual))
  if(print.results) {
    resultsVector<-c(AICValue,lnLValue,exp(par))
    names(resultsVector)<-c("AIC","lnL",msIndividualParameters(migrationIndividual))
    print(resultsVector)
  }
  return(AICValue)
}

searchContinuousModelSpace<-function(p, migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, maxParameterValue=10000, nTrees=1 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=FALSE, debug=FALSE,method="nlminb",itnmax=NULL,...) {
  modelID<-returnModel(p,migrationArrayMap)
  if(print.results) {
    resultVector<-c(modelID,p)
    names(resultVector)<-c("migrationArryIndividualID","collapseMatrix.number", "n0multiplierMap.number","migrationArray.number")
    print(resultVector)
  }
  if(is.na(modelID)) {
    return(badAIC)
  }
  else {
    paramNames<-msIndividualParameters(migrationArray[[modelID]])
    startingVals<-log(rep(1,length(paramNames))) #going to optimize in log space
    if(debug) {
      print(startingVals) 
    }
    searchResults<-optimx(par=startingVals, fn=returnAIC, method=method, migrationIndividual=migrationArray[[modelID]], migrationArrayMap=migrationArrayMap, migrationArray=migrationArray, popVector=popVector, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, print.ms.string=print.ms.string, print.results=print.results, itnmax=itnmax, debug=debug, control=list(iter.max=itnmax, maxit=itnmax,iter.lim=itnmax),...)
    return(searchResults$fvalues)
  }
}

searchContinuousModelSpaceOptim<-function(p, migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, maxParameterValue=10000, nTrees=1 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=FALSE, debug=FALSE,method="nlminb",itnmax=NULL,...) {
  modelID<-returnModel(p,migrationArrayMap)
  if(print.results) {
    resultVector<-c(modelID,p)
    names(resultVector)<-c("migrationArryIndividualID","collapseMatrix.number", "n0multiplierMap.number","migrationArray.number")
    print(resultVector)
  }
  if(is.na(modelID)) {
    return(badAIC)
  }
  else {
    paramNames<-msIndividualParameters(migrationArray[[modelID]])
    startingVals<-log(rep(1,length(paramNames))) #going to optimize in log space
    if(debug) {
      print(startingVals) 
    }
    searchResults<-optim(par=startingVals, fn=returnAIC, method=method, control=list(maxit=itnmax), migrationIndividual=migrationArray[[modelID]], migrationArrayMap=migrationArrayMap, migrationArray=migrationArray, popVector=popVector, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, print.ms.string=print.ms.string, print.results=print.results, itnmax=itnmax, debug=debug, ...)
    return(searchResults$value)
  }
}

searchDiscreteModelSpace<-function(migrationArrayMap, migrationArray, popVector, print.ms.string=FALSE, badAIC=100000000000000, maxParameterValue=10000, nTrees=1,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, debug=FALSE,  method="newuoa", itnmax=NULL,pop.size=50, print.results=FALSE, ...) {
  Domains<-matrix(ncol=2,nrow=3)
  Domains[1,]<-range(migrationArrayMap$collapseMatrix.number)
  Domains[2,]<-range(migrationArrayMap$n0multiplierMap.number)
  Domains[3,]<-range(migrationArrayMap$migrationArray.number)
  
  results<-genoud(searchContinuousModelSpace,nvars=3, max=FALSE,starting.values=c(1,1,1), MemoryMatrix=TRUE, boundary.enforcement=2, data.type.int=TRUE, Domains=Domains, migrationArrayMap=migrationArrayMap, migrationArray=migrationArray, popVector=popVector, print.ms.string=print.ms.string, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, debug=debug, method=method,itnmax=itnmax, pop.size=pop.size, print.results=print.results,...)
  return(results)
}

searchDiscreteModelSpaceOptim<-function(migrationArrayMap, migrationArray, popVector, print.ms.string=FALSE, badAIC=100000000000000, maxParameterValue=10000, nTrees=1,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, debug=FALSE,  method="newuoa", itnmax=NULL,pop.size=50, print.results=FALSE, ...) {
  Domains<-matrix(ncol=2,nrow=3)
  Domains[1,]<-range(migrationArrayMap$collapseMatrix.number)
  Domains[2,]<-range(migrationArrayMap$n0multiplierMap.number)
  Domains[3,]<-range(migrationArrayMap$migrationArray.number)
  
  results<-genoud(searchContinuousModelSpaceOptim,nvars=3, max=FALSE,starting.values=c(1,1,1), MemoryMatrix=TRUE, boundary.enforcement=2, data.type.int=TRUE, Domains=Domains, migrationArrayMap=migrationArrayMap, migrationArray=migrationArray, popVector=popVector, print.ms.string=print.ms.string, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, debug=debug, method=method,itnmax=itnmax, pop.size=pop.size, print.results=print.results,...)
  return(results)
}
