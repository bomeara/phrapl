if(getRversion() >= "2.15.1")  utils::globalVariables(c("ncores", "popVector", "maxK","migrationArray", "migrationArrayMap"))

TruncateToCells<-function(x, numCells, minVal, maxVal) {
		if(x<=numCells && min(maxVal, minVal+x)>=minVal) {
			return(seq(from=minVal, to=min(maxVal, minVal+x), by=1))
		} else {
			return(NA)	
		}
}

RowComplete<-function(x, allowableMaxInitial) {
	if(length(unique(x))==length(seq(from=min(x, allowableMaxInitial), to=max(x), by=1))) {
		return(TRUE)	
	}	
	return(FALSE)
}

AllParamCombinations<-function(numCells, minVal=0, maxVal=3, allowableMaxInitial=1) {
	#first cell can be minval...1 (so, usually, 0 or 1). Next cell can be 0, 1, 2 (as long as maxval>=2), but the second cell can be no more than one greater than all previous cells. and so forth
	#so values can be 00, 01, [but not 02], 10, 11, 12
	#each row is a possible combination
	#allowableMaxInitial allows you to have a min value of zero but let the first index be 1 if you want
	initial.grid<-expand.grid(TruncateToCells(1, numCells, minVal, maxVal), 
		TruncateToCells(2, numCells, minVal, maxVal),
		TruncateToCells(3, numCells, minVal, maxVal),
		TruncateToCells(4, numCells, minVal, maxVal),
		TruncateToCells(5, numCells, minVal, maxVal),
		TruncateToCells(6, numCells, minVal, maxVal),
		TruncateToCells(7, numCells, minVal, maxVal),
		TruncateToCells(8, numCells, minVal, maxVal),
		TruncateToCells(9, numCells, minVal, maxVal),
		TruncateToCells(10, numCells, minVal, maxVal),
		TruncateToCells(11, numCells, minVal, maxVal),
		TruncateToCells(12, numCells, minVal, maxVal),
		TruncateToCells(13, numCells, minVal, maxVal),
		TruncateToCells(14, numCells, minVal, maxVal),
		TruncateToCells(15, numCells, minVal, maxVal),
		TruncateToCells(16, numCells, minVal, maxVal),
		TruncateToCells(17, numCells, minVal, maxVal),
		TruncateToCells(18, numCells, minVal, maxVal),
		TruncateToCells(19, numCells, minVal, maxVal),
		TruncateToCells(20, numCells, minVal, maxVal),
		TruncateToCells(21, numCells, minVal, maxVal),
		TruncateToCells(22, numCells, minVal, maxVal),
		TruncateToCells(23, numCells, minVal, maxVal),
		TruncateToCells(24, numCells, minVal, maxVal),
		TruncateToCells(25, numCells, minVal, maxVal),
		TruncateToCells(26, numCells, minVal, maxVal),
		TruncateToCells(27, numCells, minVal, maxVal),
		TruncateToCells(28, numCells, minVal, maxVal),
		TruncateToCells(29, numCells, minVal, maxVal),
		TruncateToCells(30, numCells, minVal, maxVal),
		TruncateToCells(31, numCells, minVal, maxVal),
		TruncateToCells(32, numCells, minVal, maxVal),
		TruncateToCells(33, numCells, minVal, maxVal),
		TruncateToCells(34, numCells, minVal, maxVal),
		TruncateToCells(35, numCells, minVal, maxVal),
		TruncateToCells(36, numCells, minVal, maxVal),
		TruncateToCells(37, numCells, minVal, maxVal),
		TruncateToCells(38, numCells, minVal, maxVal),
		TruncateToCells(39, numCells, minVal, maxVal),
		TruncateToCells(40, numCells, minVal, maxVal),
		TruncateToCells(41, numCells, minVal, maxVal),
		TruncateToCells(42, numCells, minVal, maxVal),
		TruncateToCells(44, numCells, minVal, maxVal)
	)[,sequence(numCells)]
	if(sum(is.na(initial.grid))>0) {
		return(matrix(nrow=0, ncol=numCells))	
	}
	initial.grid<-initial.grid[which(RowMin(initial.grid)>=minVal),]
	initial.grid<-initial.grid[which(RowMax(initial.grid)<=maxVal),]
	initial.grid<-initial.grid[which(initial.grid[,1]<=allowableMaxInitial),]
	initial.grid<-initial.grid[which(apply(initial.grid, 1, RowComplete, allowableMaxInitial = allowableMaxInitial)),]
	initial.grid<-unname(as.matrix(initial.grid))
	return(initial.grid)
}

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
		if (ColCountLast(abs(popIntervalsList[[i]]$collapseMatrix),1)==0) {
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
			if (length(which(abs(lastGen)==1))>0) {
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
			rowMapping<-c(min(which(abs(lastGen)==1)),which(lastGen==0))
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

#This function will build a graph connecting populations with migration and pop coalescence. It then tests for connectivity
CheckFiniteCoalescence <- function(migrationIndividual, forceCoalescenceBasalRegime=TRUE) {
	g <- graph.empty() + vertices(sequence(dim(migrationIndividual$collapseMatrix)[1]))
	for (interval in sequence(dim(migrationIndividual$collapseMatrix)[2])) {
		merges<-which(migrationIndividual$collapseMatrix[,interval]>0) 
		if(length(merges)>1) {
			for (i in c(2:length(merges))) {
				g[merges[1], merges[i]] <- 1 #create an edge	
			}				
		}
	}
	for (i in sequence(dim(migrationIndividual$migrationArray)[1])) {
		for (j in sequence(dim(migrationIndividual$migrationArray)[2])) {
			for (k in sequence(dim(migrationIndividual$migrationArray)[3])) {
				if(!is.na(migrationIndividual$migrationArray[i, j, k])) {
					if(migrationIndividual$migrationArray[i, j, k]>0) {
						g[i, j]<-1	
					}
				}
			}	
		}
	}
	numberWithoutConnection <- sum(!is.finite(RowMax(shortest.paths(g, mode="all"))))
	finite = TRUE
	if(numberWithoutConnection>0) {
		finite=FALSE
	}
	if(finite && forceCoalescenceBasalRegime) { #here, we want to throw out models where in the rootmost regime, there is no possibility of coalescence (i.e., start with three pops with migration, two collapse, and then have no migration between the remaining 2. MS sometimes gets stuck in these cases: an allele that missed the coalescence in the first regime will just travel down the two remaining populations forever vainly seeking to coalesce with its mate in the other population, and MS doesn't detect this)
		last.regime <- dim(migrationIndividual$migrationArray)[3]
		last.populations <- which(apply(!apply(migrationIndividual$migrationArray[,, last.regime], 1, is.na), 1, sum)>0)
		g.last<- graph.empty() + vertices(sequence(length(last.populations)))
		relevant.collapseMatrix <- migrationIndividual$collapseMatrix[last.populations, last.regime]
		merges<-which(relevant.collapseMatrix>0) 
		if(length(merges)>1) {
			for (i in c(2:length(merges))) {
				g.last[merges[1], merges[i]] <- 1 #create an edge	
			}				
		}
		relevant.migrationArray <- migrationIndividual$migrationArray[last.populations, last.populations, last.regime]
		for (i in sequence(dim(relevant.migrationArray)[1])) {
			for (j in sequence(dim(relevant.migrationArray)[2])) {
				if(!is.na(relevant.migrationArray[i, j])) {
					if(relevant.migrationArray[i, j]>0) {
						g.last[i, j]<-1
					}
				}
			}
			
		}
		numberWithoutConnection <- sum(!is.finite(RowMax(shortest.paths(g.last, mode="all"))))
		if(numberWithoutConnection>0) {
			finite=FALSE
		}
	}
	return(finite)
}


GenerateIntervals<-function(popVector) {
#popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
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

# GenerateIntervalsNoCollapse<-function(popVector,maxK) {
# #popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
# #maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
  # nPop <- length(popVector)
	# firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	# firstIntervals<-firstIntervals[,ColMax(firstIntervals)==0] #we only want ones with no collapse (like an island model)
	# popIntervalsList<-list()
	# popIntervalsList[[1]]<-Popinterval(as.matrix(firstIntervals,ncol=1))
	# popIntervalsList<-CompleteIntervals(UpdateCompletes(popIntervalsList))
	# return(popIntervalsList)
# }

# GenerateIntervalsFullyResolvedCollapse<-function(popVector,maxK) {
# #popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
# #maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
  # nPop <- length(popVector)
	# firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	# firstIntervals<-firstIntervals[,ColMax(firstIntervals)<=1] #we are okay with having all zeros: no population collapse
	# firstIntervals<-firstIntervals[,ColCountIf(firstIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
	# popIntervalsList<-list()
	# for (i in 1:dim(firstIntervals)[2]) {
		# popIntervalsList[[i]]<-Popinterval(as.matrix(firstIntervals[,i],ncol=1))
	# }
	# popIntervalsList<-CompleteIntervals(UpdateCompletes(popIntervalsList))
  # intervalsToDelete<-c()
  # for (i in sequence(length(popIntervalsList))) {
    # focalCollapseMatrix<-popIntervalsList[[i]]$collapseMatrix
    # if ((min(ColMax(focalCollapseMatrix))==0) || (dim(focalCollapseMatrix)[1]!=(1+dim(focalCollapseMatrix)[2]))) {
       # intervalsToDelete<-append(intervalsToDelete,i)
    # }
  # }
  # popIntervalsList<-popIntervalsList[-intervalsToDelete]
	# return(popIntervalsList)
# }

# GenerateIntervalsSpeciesDelimitation<-function(popVector,maxK) {
# #popVector is samples from each pop: c(5,6,2) means 5 samples from pop 1, 6 from pop 2, and 2 from pop3
# #maxK is the maximum number of free parameters you want. By default, allows one free parameter for every 20 samples
  # nPop <- length(popVector)
	# firstIntervals<-blockparts(c(1:nPop),nPop,include.fewer=TRUE)
	# firstIntervals<-firstIntervals[,ColMax(firstIntervals)==1] #we are not ok with all zeros
	# firstIntervals<-firstIntervals[,ColCountIf(firstIntervals,1)!=1] #we want intervals that have at least two ones. What does it mean to have one lineage coalesce with itself?
	# firstIntervals <- cbind(firstIntervals, -firstIntervals) #where -1 = merge into one species at present
	# popIntervalsList<-list()
	# for (i in 1:dim(firstIntervals)[2]) {
		# popIntervalsList[[i]]<-Popinterval(as.matrix(firstIntervals[,i],ncol=1))
	# }
	# popIntervalsList<-CompleteIntervals(UpdateCompletes(popIntervalsList))
  # intervalsToDelete<-c()
  # for (i in sequence(length(popIntervalsList))) {
    # focalCollapseMatrix<-popIntervalsList[[i]]$collapseMatrix
    # if ((min(ColMax(focalCollapseMatrix))==0) || (dim(focalCollapseMatrix)[1]!=(1+dim(focalCollapseMatrix)[2]))) {
       # intervalsToDelete<-append(intervalsToDelete,i)
    # }
  # }
  # popIntervalsList<-popIntervalsList[-intervalsToDelete]
	# return(popIntervalsList)
# }


#the basic idea here is that at each population in each time interval there is a n0multiplier. These can all be set to the same value, allowed to vary, or assigned in clumps (i.e., pops 1, 3, and 6 have the same n0multiplier value)
#this generates all such mappings, subject to staying within the maximum number of free parameters
GenerateN0multiplierIndividuals<-function(popVector,popIntervalsList=GenerateIntervals(popVector),maxK=SetMaxK(popVector), maxN0K=Inf) {
	n0multiplierIndividualsList<-list()
	for (i in sequence(length(popIntervalsList))) {
		n0multiplierMapTemplate<-1+0*popIntervalsList[[i]]$collapseMatrix  #will have all the populations, all with either NA or 1
		numLineages=sum(n0multiplierMapTemplate,na.rm=TRUE)
		possibleMappings<-AllParamCombinations(numCells = numLineages, minVal=1, maxVal=max(1,min(maxN0K, maxK-KPopInterval(popIntervalsList[[i]]))), allowableMaxInitial=1)
		for (mappingIndex in sequence(dim(possibleMappings)[1])) {
			thisMapping<-possibleMappings[mappingIndex,]
			n0multiplierMap<-n0multiplierMapTemplate
			whichPositions <- which(n0multiplierMap==1)
			n0multiplierMap[whichPositions]<-thisMapping
			n0multiplierIndividualsList[[length(n0multiplierIndividualsList)+1]]<-N0multiplierindividual(popIntervalsList[[i]]$collapseMatrix, popIntervalsList[[i]]$complete, n0multiplierMap)
		}
	}
	return(n0multiplierIndividualsList)
}

#now we will generate all possible assignments of pairwise migration. Again, we want to keep the total number of free parameters (times, n0multipliers, migration rates) under our chosen max
#allow a model where migrations change anywhere along branch, or only at coalescent nodes? The problem with the latter is that you cannott fit some reasonable models: i.e., two populations persisting through time. Problem with the former is parameter space
GenerateMigrationIndividuals<-function(popVector, maxK=SetMaxK(popVector), maxMigrationK=2, maxN0K=1, forceSymmetricalMigration=TRUE, forceTree=FALSE, verbose=FALSE) {
	popIntervalsList<-GenerateIntervals(popVector)
	n0multiplierIndividualsList<-GenerateN0multiplierIndividuals(popVector,popIntervalsList,maxK=maxK, maxN0K=maxN0K)
	migrationIndividualsList<-list()
	for (i in sequence(length(n0multiplierIndividualsList))) {
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
							if(!forceSymmetricalMigration) {
								if (!is.na(collapseMatrix[fromPop,interval]) && !is.na(collapseMatrix[toPop,interval])) {
									migrationTemplate[fromPop,toPop,interval]<-1
								}
							} else {
								if (!is.na(collapseMatrix[fromPop,interval]) && !is.na(collapseMatrix[toPop,interval]) && toPop>fromPop) { #so only fill top diagonal
									migrationTemplate[fromPop,toPop,interval]<-1
								}								
							}
						}
					}
				}
			}
			numPairs=sum(migrationTemplate,na.rm=TRUE)
			mappingMaxKAllowed <- min(maxK - ( KCollapseMatrix(collapseMatrix) + KN0multiplierMap(n0multiplierMap) ), maxMigrationK)
			
			evaluateMapping<-TRUE
			if (mappingMaxKAllowed<=0) {
				evaluateMapping<-FALSE #since there's no way to do this and not have too many parameters
			}
			allMappings<-AllParamCombinations(numPairs, 0, mappingMaxKAllowed, 1)
			if(verbose==TRUE) {
				print(paste("  there are ", dim(allMappings)[1], " migration mappings to try", sep=""))	
			}
			thisMapping<-allMappings[1,]
			numevaluations=1
			while(evaluateMapping) {
				migrationArray<-migrationTemplate
				whichPositions <- which(migrationArray==1)
				migrationArray[whichPositions]<-thisMapping
				if(forceSymmetricalMigration) {
					for (interval in 1:numSteps) {
						for (fromPop in 1:numFinalPops) {
							for (toPop in 1:numFinalPops) {
								if (toPop<fromPop) {
									migrationArray[fromPop,toPop,interval]<-migrationArray[toPop,fromPop,interval]
								}
							}
						}
					}
				}
				newIndividual<-Migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArray)	
				if(CheckFiniteCoalescence(newIndividual) && KAll(newIndividual)<maxK) {
					migrationIndividualsList[[length(migrationIndividualsList)+1]]<-newIndividual
				}			
				numevaluations<-numevaluations+1
				if(numevaluations<=dim(allMappings)[1]) {
					thisMapping<-allMappings[numevaluations,]	
				} else {
					evaluateMapping<-FALSE	
				}
			}		
		}
	}
	if(forceTree) {
		migrationIndividualsList<-FilterForFullyResolved(migrationIndividualsList)
	}
	return(migrationIndividualsList)
}

# #this will create a set of models with no populations merging. However, not all populations need to have migration to every other population
# GenerateMigrationIndividualsNoCollapseAllowNoMigration<-function(popVector,n0multiplierIndividualsList=GenerateN0multiplierIndividuals(popVector,popIntervalsList=GenerateIntervalsNoCollapse(popVector,maxK=SetMaxK(popVector)),maxK=SetMaxK(popVector)), maxK=SetMaxK(popVector), verbose=FALSE,file=NULL) {
  # migrationArray<-GenerateMigrationIndividualsAllowNoMigration(popVector,n0multiplierIndividualsList, maxK, verbose=verbose, file=NULL)
  # migrationIndividualsToKill<-c()
  # for (i in sequence(length(migrationArray))) {
     # if (RowMax(migrationArray[[i]]$migrationArray[, , 1])==0 && RowMax(migrationArray[[i]]$migrationArray[, , 1])==0) {
        # migrationIndividualsToKill<-append( migrationIndividualsToKill,i)
     # }
  # }
  # migrationArray<-migrationArray[-1*migrationIndividualsToKill]
  # if (!is.null(file)) {
    # save(migrationArray,maxK,popVector,file=file,compress=TRUE)
  # }
  # return(migrationArray)
# }

# #this will create a set of models where populations are linked only by a bifurcating tree
# GenerateMigrationIndividualsFullyResolvedCollapseAllowNoMigration<-function(popVector,n0multiplierIndividualsList=GenerateN0multiplierIndividuals(popVector,popIntervalsList=GenerateIntervalsFullyResolvedCollapse(popVector,maxK=SetMaxK(popVector)),maxK=SetMaxK(popVector)), maxK=SetMaxK(popVector), verbose=FALSE,file=NULL) {
  # migrationArray<-GenerateMigrationIndividualsAllowNoMigration(popVector,n0multiplierIndividualsList, maxK, verbose=verbose, file=NULL)
  # migrationIndividualsToKill<-c()
  # for (i in sequence(length(migrationArray))) {
     # if (RowMax(migrationArray[[i]]$migrationArray[, , 1])==0 && RowMax(migrationArray[[i]]$migrationArray[, , 1])==0) {
        # migrationIndividualsToKill<-append( migrationIndividualsToKill,i)
     # }
  # }
  # migrationArray<-migrationArray[-1*migrationIndividualsToKill]
  # if (!is.null(file)) {
    # save(migrationArray,maxK,popVector,file=file,compress=TRUE)
  # }
  # return(migrationArray)
# }


FilterForFullyResolved <- function(migrationArray) {
	migrationIndividualsToKill<-c()
   for (i in sequence(length(migrationArray))) {
      if(min(ColMax(migrationArray[[i]]$collapseMatrix))==0 || dim(migrationArray[[i]]$collapseMatrix)[2]!=(dim(migrationArray[[i]]$collapseMatrix)[1]-1)) {
         migrationIndividualsToKill<-append( migrationIndividualsToKill,i)
      }
   }
   migrationArray<-migrationArray[-1*migrationIndividualsToKill]
   return(migrationArray)
}

# #this is like generateMigrationIndividuals, but allows some migration rates to be set to 0 with no penalty in terms of number of free parameters
# GenerateMigrationIndividualsAllowNoMigration<-function(popVector,n0multiplierIndividualsList=GenerateN0multiplierIndividuals(popVector,popIntervalsList=GenerateIntervals(popVector,maxK),maxK), maxK, verbose=FALSE, file=NULL) {
	# migrationIndividualsList<-list()
	# for (i in 1:length(n0multiplierIndividualsList)) {
		# if (verbose==TRUE) {
			# print(paste("doing ",i,"/",length(n0multiplierIndividualsList), " (migration individuals so far = ",length(migrationIndividualsList),")"))
		# }
		# collapseMatrix<-n0multiplierIndividualsList[[i]]$collapseMatrix
		# n0multiplierMap<-n0multiplierIndividualsList[[i]]$n0multiplierMap
		# numFinalPops<-dim(collapseMatrix)[1]
		# numSteps<-dim(collapseMatrix)[2]
		# if ((KCollapseMatrix(collapseMatrix) + KN0multiplierMap(n0multiplierMap) )<=maxK) {
			# migrationTemplate<-array(data=NA,dim=c(numFinalPops,numFinalPops,numSteps),dimnames=c("from","to","generation"))
			# for (interval in 1:numSteps) {
				# for (fromPop in 1:numFinalPops) {
					# for (toPop in 1:numFinalPops) {
						# if (fromPop!=toPop) {
							# if (!is.na(collapseMatrix[fromPop,interval]) && !is.na(collapseMatrix[toPop,interval])) {
								# migrationTemplate[fromPop,toPop,interval]<-1
							# }
						# }
					# }
				# }
			# }
			# numPairs=sum(migrationTemplate,na.rm=TRUE)
			# mappingMaxKAllowed <- maxK - ( KCollapseMatrix(collapseMatrix) + KN0multiplierMap(n0multiplierMap) )
			# evaluateMapping<-TRUE
			# if (mappingMaxKAllowed<0) {
				# evaluateMapping<-FALSE #since there's no way to do this and not have too many parameters
			# }
			# thisMapping<-firstcomposition(numPairs)		
			# while(evaluateMapping) {
				# thisMapping<-c(thisMapping,rep(0,numPairs-length(thisMapping))) #adds trailing zeros
				# thisMappingOrig<-thisMapping
				# if ((length(which(thisMapping>0)) - 1 ) <=mappingMaxKAllowed) { #note the -1: we take away one free param because we allow zero migration rates
					# migrationArray<-migrationTemplate
					# whichPositions <- which(migrationArray==1)
					# for (positionIndex in 1:length(whichPositions)) {
						# position=whichPositions[positionIndex]
						# paramPosition<-which(thisMapping>0)[1]
						# migrationArray[position]=paramPosition #the position of the first parameter
						# thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we have used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
					# }
					# if (is.na(max(migrationArray,na.rm=TRUE))) {
						# migrationIndividualsList[[length(migrationIndividualsList)+1]]<-Migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArray)				
						# #print(paste("Just created migration individual ",length(migrationIndividualsList)))
					# }
					# else {
						# for (paramToMakeZero in sequence(max(migrationArray,na.rm=TRUE))) {
							# migrationArrayModified<-migrationArray
							# migrationArrayModified[which(migrationArrayModified==paramToMakeZero)]<-0
							# migrationArrayModified[which(migrationArrayModified>paramToMakeZero)]<-migrationArrayModified[which(migrationArrayModified>paramToMakeZero)]-1
							# migrationIndividualsList[[length(migrationIndividualsList)+1]]<-Migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArrayModified)
							# #print(paste("Just created migration individual ",length(migrationIndividualsList)))
						# }
					# }
				# }
				# if (!islastcomposition(thisMappingOrig,restricted=FALSE)) {
					# thisMapping<-nextcomposition(thisMappingOrig,restricted=FALSE)
				# }
				# else {
					# evaluateMapping<-FALSE
				# }
			# }
		# }
	# }
  # if (!is.null(file)) {
    # migrationArray<-migrationIndividualsList
    # save(migrationArray,maxK,popVector,file=file,compress=TRUE)
  # }
	# return(migrationIndividualsList)
# }

# GenerateExpansionMultiplierIndividuals<-function(popVector,migrationIndividualsList=GenerateMigrationIndividualsAllowNoMigration(popVector,n0multiplierIndividualsList=GenerateN0multiplierIndividuals(popVector,popIntervalsList=GenerateIntervals(popVector,maxK),maxK),maxK),maxK) {
	# expansionmultiplierIndividualsList<-list()
	# for (i in 1:length(migrationIndividualsList)) {
		# expansionmultiplierMapTemplate<-1+0*migrationIndividualsList[[i]]$collapseMatrix  #will have all the populations, all with either NA or 1
		# numLineages=sum(expansionmultiplierMapTemplate,na.rm=TRUE)
		# possibleMappings<-compositions(numLineages)
		# for (mappingIndex in 1:dim(possibleMappings)[2]) {
			# thisMapping<-possibleMappings[,mappingIndex]
			# if ((length(which(thisMapping>0))+KPopInterval(migrationIndividualsList[[i]]) )<=maxK) { #only do it on those mappings that have not too many free parameters
				# expansionmultiplierMap<-expansionmultiplierMapTemplate
				# whichPositions <- which(expansionmultiplierMap==1)
				# for (positionIndex in 1:length(whichPositions)) {
					# position=whichPositions[positionIndex]
					# paramPosition<-which(thisMapping>0)[1]
					# expansionmultiplierMap[position]=paramPosition #the position of the first parameter
					# thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we have used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
				# }
				# expansionmultiplierIndividualsList[[length(expansionmultiplierIndividualsList)+1]]<-Expansionmultiplierindividual(migrationIndividualsList[[i]]$collapseMatrix, migrationIndividualsList[[i]]$complete, expansionmultiplierMap)
			# }
		# }
	# }
	# return(expansionmultiplierIndividualsList)
# }


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

#if "popAssignments" has not been specified, then iteratively randomly sample "nIndividualsDesired" from the 
#entire dataset "attemptsCutoff" times until at least "minPerPop" individuals are represented per population 
#within "toRetain"	
TaxaToRetain<-function(assignFrame,nIndividualsDesired,minPerPop=1,attemptsCutoff=100000,popAssignments=NULL) {
	if(is.null(popAssignments)) {
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
		#if popAssignments has been specified, then randomly sample according to the population-specific sample sizes
		#indicated within popAssignments
		nIndividualsDesired<-sum(popAssignments)
		toRetain <- array()
		#for each population...
		for(numLevel in 1:length(unique(assignFrame[,1])))
		{
			thisPopSamples <- assignFrame[which(assignFrame[,1] == unique(assignFrame[,1])[numLevel]),][[2]]
			#...if the there is more than one individual in the population...
			if (length(thisPopSamples) > 1)			{
				#...randomly sample "popAssignments[numLevel]" individuals and add them to "toRetain"
				toRetainTemp <- sample(thisPopSamples,popAssignments[numLevel],replace=FALSE)
				toRetain <- c(toRetain,toRetainTemp)[!is.na(c(toRetain,toRetainTemp))]
			} else {
				#But if there is only one individual in the population (and one sample is to be drawn)...
				if(popAssignments[numLevel] == 1) {
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
	taxaToDrop<-allTaxa[-as.numeric(taxaRetained)]
	return(taxaToDrop)
}

PrepSubsampling<-function(assignmentsGlobal, observedTrees,popAssignments,subsamplesPerGene,nIndividualsDesired=NULL,minPerPop=1,
outgroup=TRUE,outgroupPrune=TRUE){
	if(is.null(nIndividualsDesired)){
		nIndividualsDesired<-sum(popAssignments[[1]])
	}
	#Create list for storing subsample trees
	phyList<-list()
	#If outgroup present, add this to the first popAssignments vector
	if(outgroup==TRUE){
		popAssignments[[1]]<-c(popAssignments[[1]],1)
		nIndividualsDesired<-nIndividualsDesired + 1
	}
	assignFrameOriginal<-assignmentsGlobal
	assignFrameOriginal<-cbind(assignFrameOriginal,c(1:nrow(assignFrameOriginal)))
	colnames(assignFrameOriginal)<-c("indiv","popLabel","indivTotal")
	phyOriginal <- observedTrees
	if(class(phyOriginal)!="multiPhylo") { #if there is only one locus...
		phyOriginal<-c(phyOriginal) #convert phy to a multiphylo class
	}
	#Create list for storing subsampled trees
	counters<-c()
	for(i in sequence(length(popAssignments))){
		phyList[[i]]<-rmtree(N=length(phyOriginal) * subsamplesPerGene,n=3)
		counters<-append(counters,1)
	}
	
	for (tree in 1:length(phyOriginal)){ #for each locus (i.e., tree)
		phy<-phyOriginal #renew phy for each locus
		assignFrame<-assignFrameOriginal	#renew assignFrame for each locus	
		tips.vec <- phy[[tree]]$tip.label #make an array of the included individuals
		phy[[tree]]$tip.label <- as.character(c(1:length(phy[[tree]]$tip.label))) #change tip labels to consecutive numbers			
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
			} else { #unless there was no matching sample in the assingment file, in which case, drop this tip from the tree
				phy[[tree]]<-drop.tip(phy[[tree]],as.character(match)) 
				cat("Warning: Tree number ",tree,"contains tip names not included in the inputted assignment file.",
					"These tips will not be subsampled.\n",file=paste(subsamplePath,"WARNINGS.txt",sep=""),append=TRUE)
			}
		}
		assignFrame<-rbind(data.frame(match.vec$popLabel,c(1:length(phy[[tree]]$tip.label)))) #convert this locus-specific one into assignFrame
		colnames(assignFrame) <- c("popLabel","indivTotal")	
		assignFrame <- assignFrame[order(assignFrame$popLabel),] #order new assignFrame by population
	
		#Re-name the tree tips so they match the assignment file order (tips belonging to population A are numbered first, B second, C third, etc)
		for(changetips in 1:length(phy[[tree]]$tip.label))
		{
			phy[[tree]]$tip.label[assignFrame$indivTotal[changetips]] <- as.numeric(changetips)
		}

		#Make new assignFrame with the new tip labels listed in order
		assignFrame <- data.frame(as.factor(assignFrame[,1]),c(1:length(phy[[tree]]$tip.label)))
		colnames(assignFrame) <- c("popLabel","indivTotal")	
		
		#Begin the subsampling
		retainedTaxaMatrix<-matrix(NA,nrow=subsamplesPerGene,ncol=nIndividualsDesired) #matrix storing subsamples from each iteration
		newphyVector<-list()
		for (rep in 1:subsamplesPerGene) {
			keepTaxa<-TaxaToRetain(assignFrame,nIndividualsDesired,minPerPop,attemptsCutoff=100000,popAssignments[[1]]) #subsample assignFrame
			retainedTaxaMatrix[rep,]<-keepTaxa
			prunedAF<-PrunedAssignFrame(assignFrame,keepTaxa)
			delTaxa<-TaxaToDrop(assignFrame,keepTaxa)
			newphy<-drop.tip(phy[[tree]],as.character(delTaxa)) #toss non-sampled individuals from the tree
			for (tipIndex in 1:length(newphy$tip.label)) { #rename tips in tree to be consecutive
				old.label<-newphy$tip.label[tipIndex]
				new.label<-as.character(which(prunedAF[,2]==old.label))
				newphy$tip.label[tipIndex]<-new.label
			}
			#If tossing outgroup
			if(outgroupPrune==TRUE){
				notFound=FALSE
				tipIndex=0
				
				#Toss out that tip which has the highest tip label
				while(notFound==FALSE){
					tipIndex=tipIndex + 1
					if(newphy$tip.label[tipIndex]==length(keepTaxa)){
						newphy<-drop.tip(newphy,tipIndex)
						notFound=TRUE
					}
				}
			}		
			prunedAF[,2]<-as.character(c(1:length(prunedAF[,2]))) #rename tips in assignFrame to be consecutive
			if(outgroupPrune==TRUE){
				prunedAF<-prunedAF[-length(prunedAF[,2]),] #prune outgroup from assignFrame
			}
			
			#Save current subsample
			phyList[[1]][[counters[1]]]<-newphy
			counters[1]<-counters[1] + 1

#			write.tree(newphy,file=paste(subsamplePath,"observed1.tre",sep=""),append=TRUE)
			#Save current subsample in master list
			newphyVector[[length(newphyVector) + 1]]<-newphy
		}
		
		#Subsample further if specified by popAssignments, printing to the tree file each time
		if(length(popAssignments) > 1){
			assignFrame<-prunedAF
			for(i in 2:(length(popAssignments))){
				currentSampleSize=sum(popAssignments[[i]])
				retainedTaxaMatrix<-matrix(NA,nrow=subsamplesPerGene,ncol=sum(popAssignments[[i]])) #matrix storing subsamples from each iteration	
				for(j in 1:subsamplesPerGene){
					keepTaxa<-TaxaToRetain(assignFrame,nIndividualsDesired=currentSampleSize,
						minPerPop,attemptsCutoff=100000,popAssignments=popAssignments[[i]]) #subsample assignFrame
					retainedTaxaMatrix[j,]<-keepTaxa
					prunedAF<-PrunedAssignFrame(assignFrame,keepTaxa)
					delTaxa<-TaxaToDrop(assignFrame,taxaRetained=keepTaxa)
					newphy<-drop.tip(newphyVector[[j]],as.character(delTaxa)) #toss non-sampled individuals from the tree
					prunedAF<-cbind(prunedAF[order(prunedAF[,1],prunedAF[,2]),],new.label=c(1:nrow(prunedAF)))
					for (tipIndex in 1:length(newphy$tip.label)) { #rename tips in tree to be consecutive
						old.label<-newphy$tip.label[tipIndex]
						new.label<-as.character(prunedAF[,3][which(prunedAF[,2]==old.label)])
						newphy$tip.label[tipIndex]<-new.label
					}
					#Save current subsample
					phyList[[i]][[counters[i]]]<-newphy
					counters[i]<-counters[i] + 1
					
#					write.tree(newphy,file=paste(subsamplePath,"observed",i,".tre",sep=""),append=TRUE)
					newphyVector[[j]]<-newphy
				}
				assignFrame<-CreateAssignment.df(popAssignments[[i]])
			}
		}
	}
	return(phyList)
}


#This function inputs 1) an assignment file that includes all samples pooled from across loci and 2) a tree file 
#containing a tree for each locus. For each locus, subsampling can be done either by iteratively sampling 
#nIndividualsDesired from the entire dataset (with a minimum sample per population specified by minPerPop), 
#or, if popAssignments is specified, can be done by sampling a specified number of individuals per population.
#A single outgroup can also be included in each subsample, and then pruned from the tree. Input and output is 
#placed in a specified subsamplePath and includes 1) a subsamplesPerGene number of tree files with subsampled trees 
#from each locus and 2) a single assignment file (if popAssignments is specified) or an assignment file for each locus 
#and replicate (if popAssignments=NULL). If more than one subsampling vector is included in popAssignments, subsampling
#is done for each subsampling size class.  
PrepSubsamplingFiles<-function(subsamplePath="./",assignFile="cladeAssignments.txt",treesFile="trees.tre",
outputFile="observed.tre",popAssignments,subsamplesPerGene,nIndividualsDesired=NULL,minPerPop=1,
outgroup=TRUE,outgroupPrune=TRUE){
	if(is.null(nIndividualsDesired)){
		nIndividualsDesired<-sum(popAssignments[[1]])
	}
	#Remove old files
	if(file.exists(paste(subsamplePath,outputFile,sep=""))){
		unlink(paste(subsamplePath,outputFile,sep=""))
	}
	#Create list for storing subsample trees
	phyList<-list()
	#If outgroup present, add this to the first popAssignments vector
	if(outgroup==TRUE){
		popAssignments[[1]]<-c(popAssignments[[1]],1)
		nIndividualsDesired<-nIndividualsDesired + 1
	}
	#Read in and prep the assignment file constructed from all loci
	assignFrameOriginal<-read.table(paste(subsamplePath,assignFile,sep=""),header=TRUE) #read in the assignemt file including all individuals
	assignFrameOriginal<-cbind(assignFrameOriginal,c(1:nrow(assignFrameOriginal)))
	colnames(assignFrameOriginal)<-c("indiv","popLabel","indivTotal")
	#Read in tree file
	phy.file <- paste(subsamplePath,treesFile,sep="")
	phyOriginal <- read.tree(file=phy.file) #read in tree from file
	if(class(phyOriginal)!="multiPhylo") { #if there is only one locus...
		phyOriginal<-c(phyOriginal) #convert phy to a multiphylo class
	}
	#Creat list for storing subsampled trees
	counters<-c()
	for(i in sequence(length(popAssignments))){
		phyList[[i]]<-rmtree(N=length(phyOriginal) * subsamplesPerGene,n=3)
		counters<-append(counters,1)
	}
	
	for (tree in 1:length(phyOriginal)){ #for each locus (i.e., tree)
		phy<-phyOriginal #renew phy for each locus
		assignFrame<-assignFrameOriginal	#renew assignFrame for each locus	
		tips.vec <- phy[[tree]]$tip.label #make an array of the included individuals
		phy[[tree]]$tip.label <- as.character(c(1:length(phy[[tree]]$tip.label))) #change tip labels to consecutive numbers			
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
			} else { #unless there was no matching sample in the assingment file, in which case, drop this tip from the tree
				phy[[tree]]<-drop.tip(phy[[tree]],as.character(match)) 
				cat("Warning: Tree number ",tree,"contains tip names not included in the inputted assignment file.",
					"These tips will not be subsampled.\n",file=paste(subsamplePath,"WARNINGS.txt",sep=""),append=TRUE)
			}
		}
		assignFrame<-rbind(data.frame(match.vec$popLabel,c(1:length(phy[[tree]]$tip.label)))) #convert this locus-specific one into assignFrame
		colnames(assignFrame) <- c("popLabel","indivTotal")	
		assignFrame <- assignFrame[order(assignFrame$popLabel),] #order new assignFrame by population
	
		#Re-name the tree tips so they match the assignment file order (tips belonging to population A are numbered first, B second, C third, etc)
		for(changetips in 1:length(phy[[tree]]$tip.label))
		{
			phy[[tree]]$tip.label[assignFrame$indivTotal[changetips]] <- as.numeric(changetips)
		}

		#Make new assignFrame with the new tip labels listed in order
		assignFrame <- data.frame(as.factor(assignFrame[,1]),c(1:length(phy[[tree]]$tip.label)))
		colnames(assignFrame) <- c("popLabel","indivTotal")	
		
		#Begin the subsampling
		retainedTaxaMatrix<-matrix(NA,nrow=subsamplesPerGene,ncol=nIndividualsDesired) #matrix storing subsamples from each iteration
		newphyVector<-list()
		for (rep in 1:subsamplesPerGene) {
			keepTaxa<-TaxaToRetain(assignFrame,nIndividualsDesired,minPerPop,attemptsCutoff=100000,popAssignments[[1]]) #subsample assignFrame
			retainedTaxaMatrix[rep,]<-keepTaxa
			prunedAF<-PrunedAssignFrame(assignFrame,keepTaxa)
			delTaxa<-TaxaToDrop(assignFrame,keepTaxa)
			newphy<-drop.tip(phy[[tree]],as.character(delTaxa)) #toss non-sampled individuals from the tree
			for (tipIndex in 1:length(newphy$tip.label)) { #rename tips in tree to be consecutive
				old.label<-newphy$tip.label[tipIndex]
				new.label<-as.character(which(prunedAF[,2]==old.label))
				newphy$tip.label[tipIndex]<-new.label
			}
			#If tossing outgroup
			if(outgroupPrune==TRUE){
				notFound=FALSE
				tipIndex=0
				
				#Toss out that tip which has the highest tip label
				while(notFound==FALSE){
					tipIndex=tipIndex + 1
					if(newphy$tip.label[tipIndex]==length(keepTaxa)){
						newphy<-drop.tip(newphy,tipIndex)
						notFound=TRUE
					}
				}
			}		
			prunedAF[,2]<-as.character(c(1:length(prunedAF[,2]))) #rename tips in assignFrame to be consecutive
			if(outgroupPrune==TRUE){
				prunedAF<-prunedAF[-length(prunedAF[,2]),] #prune outgroup from assignFrame
			}
			
			#Save current subsample
			phyList[[1]][[counters[1]]]<-newphy
			counters[1]<-counters[1] + 1

#			write.tree(newphy,file=paste(subsamplePath,"observed1.tre",sep=""),append=TRUE)
			#Save current subsample in master list
			newphyVector[[length(newphyVector) + 1]]<-newphy
			
			if(is.null(popAssignments)){
				if (!file.exists(paste(subsamplePath,"assignments",sep=""))){
					dir.create(paste(subsamplePath,"assignments",sep=""))
				}
				write.table(prunedAF,file=paste(subsamplePath,"assignments/locus",tree,".assign",rep,".txt",sep=""),quote=FALSE,sep="\t",
				row.names=FALSE,col.names=FALSE,append=FALSE)
			}
		}
		
		#Subsample further if specified by popAssignments, printing to the tree file each time
		if(length(popAssignments) > 1){
			assignFrame<-prunedAF
			for(i in 2:(length(popAssignments))){
				currentSampleSize=sum(popAssignments[[i]])
				retainedTaxaMatrix<-matrix(NA,nrow=subsamplesPerGene,ncol=sum(popAssignments[[i]])) #matrix storing subsamples from each iteration	
				for(j in 1:subsamplesPerGene){
					keepTaxa<-TaxaToRetain(assignFrame,nIndividualsDesired=currentSampleSize,
						minPerPop,attemptsCutoff=100000,popAssignments=popAssignments[[i]]) #subsample assignFrame
					retainedTaxaMatrix[j,]<-keepTaxa
					prunedAF<-PrunedAssignFrame(assignFrame,keepTaxa)
					delTaxa<-TaxaToDrop(assignFrame,taxaRetained=keepTaxa)
					newphy<-drop.tip(newphyVector[[j]],as.character(delTaxa)) #toss non-sampled individuals from the tree
					prunedAF<-cbind(prunedAF[order(prunedAF[,1],prunedAF[,2]),],new.label=c(1:nrow(prunedAF)))
					for (tipIndex in 1:length(newphy$tip.label)) { #rename tips in tree to be consecutive
						old.label<-newphy$tip.label[tipIndex]
						new.label<-as.character(prunedAF[,3][which(prunedAF[,2]==old.label)])
						newphy$tip.label[tipIndex]<-new.label
					}
					#Save current subsample
					phyList[[i]][[counters[i]]]<-newphy
					counters[i]<-counters[i] + 1
					
#					write.tree(newphy,file=paste(subsamplePath,"observed",i,".tre",sep=""),append=TRUE)
					newphyVector[[j]]<-newphy
				}
				assignFrame<-CreateAssignment.df(popAssignments[[i]])
			}
		}
	}
	
	#Print subsampled trees
	for(m in 1:length(phyList)){
		write.tree(phyList[[m]],file=paste(subsamplePath,outputFile,sep=""),append=TRUE)
	}
	
	#Print assign frames
#	for(k in 1:length(popAssignments)){
#		if(outgroup==TRUE){
#			popAssignments[[1]]<-popAssignments[[1]][-length(popAssignments[[1]])]
#		}
#		write.table(CreateAssignment.df(popAssignments[[k]]),file=paste(subsamplePath,"assign",k,".txt",sep=""),
#			quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=FALSE)
#	}
#	unlink(paste(subsamplePath,"cladeAssignments.txt",sep=""))
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
	return(assignFrame[as.numeric(taxaRetained),])
}

#Get weights for each subsample based on the number of matches per permutation
GetPermutationWeightsAcrossSubsamples<-function(popAssignments,observedTrees){
	subsampleWeights<-data.frame(weight=matrix(1.0, ncol=1, nrow=length(observedTrees)))
	subsampleSizeCounter<-1
	for(y in sequence(length(observedTrees))){
		subsampleWeights[y,1]<-GetPermutationWeights(phy=observedTrees[[y]],popVector=popAssignments[[subsampleSizeCounter]])
		cat(paste("Finished with ",y," out of ",length(observedTrees)," trees\n",sep=""))
		if(y%%(length(observedTrees) / length(popAssignments)) == 0){
			subsampleSizeCounter<-subsampleSizeCounter + 1
		}
	}
	return(subsampleWeights)
}

#Get weights for a given tree based on the number of matches per permutation
GetPermutationWeights<-function(phy,popVector){
	assignFrame<-CreateAssignment.df(popVector)
	popLabels<-unique(assignFrame[,1])
	#Make list of label permutations for each population
	tips.permute<-list()
	try.label <- c()
	for(i in 1:length(popLabels)){
		tips.permute[[length(tips.permute) + 1]]<-as.matrix(perms(length(assignFrame[,2][which(assignFrame[,1] == 
			popLabels[i])])))
		tips.permute[[i]]<-tips.permute[[i]] + (((assignFrame[,2][which(assignFrame[,1] == popLabels[i])])[1]) - 1)
	}
	numberOfPermutations<-0
	phyOriginal<-phy
	cladesMS<-GetCladesQuickly(phy)
	matches<-0
	if(length(popLabels) > 10){
		stop("This function can only handle 10 populations")
	}

	for(j in 1:ncol(tips.permute[[1]])){
		try.label<-tips.permute[[1]][,j]
		#Define index of tips in this population
		pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[1])]))
		#Order them in the same way each time
		phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
		#Rename according to the current permutation
		phy$tip.label[pop.tips]<-try.label
					
		if(length(popLabels) > 1){	
			for(k in 1:ncol(tips.permute[[2]])){
				try.label<-tips.permute[[2]][,k]
				pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[2])]))
				phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
				phy$tip.label[pop.tips]<-try.label				

				if(length(popLabels) > 2){	
					for(l in 1:ncol(tips.permute[[3]])){
						try.label<-tips.permute[[3]][,l]
						pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[3])]))
						phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
						phy$tip.label[pop.tips]<-try.label	

						if(length(popLabels) > 3){
							for(m in 1:ncol(tips.permute[[4]])){
								try.label<-tips.permute[[4]][,m]
								pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[4])]))
								phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
								phy$tip.label[pop.tips]<-try.label	
							
								if(length(popLabels) > 4){	
									for(n in 1:ncol(tips.permute[[5]])){
										try.label<-tips.permute[[5]][,n]
										pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[5])]))
										phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
										phy$tip.label[pop.tips]<-try.label	
									
										if(length(popLabels) > 5){	
											for(o in 1:ncol(tips.permute[[6]])){
												try.label<-tips.permute[[6]][,o]
												pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[6])]))
												phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
												phy$tip.label[pop.tips]<-try.label	
											
												if(length(popLabels) > 6){	
													for(p in 1:ncol(tips.permute[[7]])){
														try.label<-tips.permute[[7]][,p]
														pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[7])]))
														phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
														phy$tip.label[pop.tips]<-try.label	
																	
														if(length(popLabels) > 7){	
															for(q in 1:ncol(tips.permute[[8]])){
																try.label<-tips.permute[[8]][,q]
																pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[8])]))
																phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
																phy$tip.label[pop.tips]<-try.label																	
														
																if(length(popLabels) > 8){
																	for(r in 1:ncol(tips.permute[[9]])){
																		try.label<-tips.permute[[9]][,r]
																		pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[9])]))
																		phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
																		phy$tip.label[pop.tips]<-try.label																			
													
																		if(length(popLabels) > 9){	
																			for(s in 1:ncol(tips.permute[[10]])){
																				try.label<-tips.permute[[10]][,s]
																				pop.tips<-which(phy$tip.label %in% as.character(assignFrame[,2][which(assignFrame[,1] == popLabels[10])]))
																				phy$tip.label[pop.tips]<-phy$tip.label[pop.tips][order(as.numeric(phy$tip.label[pop.tips]))]
																				phy$tip.label[pop.tips]<-try.label				
																			
																				if(length(popLabels) == 10){
																					cladesGene<-GetCladesQuickly(phy)
																					matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
																					numberOfPermutations<-numberOfPermutations + 1
																				}
																			}
																		}
																		if(length(popLabels) == 9){
																			cladesGene<-GetCladesQuickly(phy)
																			matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
																			numberOfPermutations<-numberOfPermutations + 1
																		}
																	}
																}
																if(length(popLabels) == 8){
																	cladesGene<-GetCladesQuickly(phy)
																	matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
																	numberOfPermutations<-numberOfPermutations + 1
																}
															}
														}
														if(length(popLabels) == 7){
															cladesGene<-GetCladesQuickly(phy)
															matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
															numberOfPermutations<-numberOfPermutations + 1
														}
													}
												}
												if(length(popLabels) == 6){
													cladesGene<-GetCladesQuickly(phy)
													matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
													numberOfPermutations<-numberOfPermutations + 1
												}
											}
										}
										if(length(popLabels) == 5){
											cladesGene<-GetCladesQuickly(phy)
											matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
											numberOfPermutations<-numberOfPermutations + 1
										}
									}
								}
								if(length(popLabels) == 4){
									cladesGene<-GetCladesQuickly(phy)
									matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
									numberOfPermutations<-numberOfPermutations + 1
								}
							}
						}
						if(length(popLabels) == 3){
							cladesGene<-GetCladesQuickly(phy)
							matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
							numberOfPermutations<-numberOfPermutations + 1
						}
					}
				}
				if(length(popLabels) == 2){
					cladesGene<-GetCladesQuickly(phy)
					matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
					numberOfPermutations<-numberOfPermutations + 1
				}
			}
		}
		if(length(popLabels) == 1){
			cladesGene<-GetCladesQuickly(phy)
			matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
			numberOfPermutations<-numberOfPermutations + 1
		}
	}
	weight<-matches / numberOfPermutations

	return(weight)
}						

##Much faster way to get a vector of clades in a tree than "GetClades (doesn't use ape's slow 
##"subtrees" function).This is currently being used with GetPermutationWeights
GetCladesQuickly <- function(phy, do.reorder=TRUE) {
	if(do.reorder) {
		phy<-reorder.phylo(phy, order="postorder")
	}
	phy$edge[which(phy$edge[,2]<=Ntip(phy)),2]<-phy$tip.label[phy$edge[which(phy$edge[,2]<=Ntip(phy)),2]]
	internal.ordered.nodes<-unique(phy$edge[,1])
	for (node.index in sequence(length(internal.ordered.nodes))) {
		phy$edge[which(phy$edge[,2]==internal.ordered.nodes[node.index]),2]<-paste(phy$edge[which(phy$edge[,1]==internal.ordered.nodes[node.index]),2], collapse="_")
	}
	clades<-sapply(c(unique(phy$edge[which(grepl("_", phy$edge[,2])),2]),paste(phy$tip.label, collapse="_")), SplitAndSort, USE.NAMES=FALSE)
	return(clades)
}

#Used with GetCladesQuickly
SplitAndSort <- function(x) {
	return(paste(sort(strsplit(x, "_")[[1]]), collapse="_"))
}


#Because these permutations take a long time, this function attempts to speed up the process
#of caluculating tree weights by randomly sampling from the possible permutations until a particular 
#level of confidence is obtained in the estimated proportion of matches. Because the number of matching 
#permutations is generally low however, this function did not end up speeding up the process. 
#On the contrary...
GetPermutationWeightsBasedOnSampling<-function(phy,popVector, desiredConfidence = 0.05, minSamples = 100){
	assignFrame<-CreateAssignment.df(popVector)
	popLabels<-unique(assignFrame[,1])
	
	#Make list of label permutations for each population
	tips.permute<-list()
	tips.permutation.order <- list()
	total.number.permutations <- 1
	seen.samples <- c()
	for(i in 1:length(popLabels)){
		tips.permute[[length(tips.permute) + 1]]<-as.matrix(perms(length(assignFrame[,2][which(assignFrame[,1] == 
			popLabels[i])])))
		tips.permute[[i]]<-tips.permute[[i]] + (((assignFrame[,2][which(assignFrame[,1] == popLabels[i])])[1]) - 1)
		tips.permutation.order[[i]] <- sample.int(dim(tips.permute[[i]])[2])
		total.number.permutations <- total.number.permutations * dim(tips.permute[[i]])[2]
	}

	numberOfPermutations<-0
	phyOriginal<-phy
	cladesMS<-GetCladesQuickly(phy)
	matches<-0
#	minRequired <- minSamples
#	while((numberOfPermutations < minSamples) || (numberOfPermutations < minRequired)) {
	while((numberOfPermutations < minSamples) || (matches < (desiredConfidence * minSamples))) {
		if(numberOfPermutations >= total.number.permutations) {
			break()
		}
		#Randomly sample permutation columns for each population (try again if that combination has been tried)
		sample.try <- rep(NA, length(tips.permute))
		while(is.na(sample.try[1])) { #we are trying to sample without replacement
			sample.try <- rep(NA, length(tips.permute))
			for (population in sequence(length(tips.permute))) {
				sample.try[population] <- sample.int(dim(tips.permute[[population]])[2],1)	
			}
			if(paste(sample.try,collapse="_") %in% seen.samples) {
				sample.try <- rep(NA, length(tips.permute))
			}else{
			}
		}
		seen.samples<-rbind(seen.samples,paste(sample.try,collapse="_"))
		
		#Create new set of tip labels
		phy <- phyOriginal
		for (population in sequence(length(tips.permute))) {
			try.label<-tips.permute[[population]][,sample.try[population]]
			phy$tip.label[order(as.numeric(phy$tip.label))][which(assignFrame[,1] == popLabels[population])] <- try.label
		}
		cat(paste(paste(as.character(phy$tip.label),collapse="_"),"\n",sep=""))	
		cladesGene<-GetCladesQuickly(phy)
		matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
		numberOfPermutations<-numberOfPermutations + 1
		proportion<-matches/numberOfPermutations
		
		#The original idea was to pose a sampling precision (e.g., 0.01) and then keep sampling until
		#the error around the estimated proportion was decreased to this point
#		minRequired<-((desiredSamplePrecision)^2) * proportion / (1-proportion)
#		minRequired<-((desiredSamplePrecision)^2) / (proportion * (1 - proportion)) 		
		weight<-matches / numberOfPermutations
	}
	
	return(weight)
}

#uses Kish's effective sample size formula
GetKishESS<-function(popVector, nsamples) {
	group.weights<-nsamples/popVector
	actual.weights<-c()
	for (i in sequence(length(group.weights))) {
		actual.weights<-c(actual.weights, rep(group.weights[i], popVector[i]))	
	}
	return(((sum(actual.weights))^2)/sum(actual.weights ^2))
}

#Uses ln likelihood
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

#Function for running seqConverter to convert nexus, fasta, or se-al files to phylip
#inputFormat can be "nex", "fasta", or "seal" (and files must use these as the file suffix)
RunSeqConverter <- function(seqConvPath=system.file("extdata","seqConverter.pl",package="phrapl"),
	inFilePath="./",inputFormat="nex"){
	filesList <- list.files(inFilePath,pattern=paste("*.",inputFormat,sep=""))
	for(numLoci in 1:length(filesList)){
		systemCall1 <-system(paste(seqConvPath," -d",paste(inFilePath,filesList[numLoci],sep="")," -ope",sep=""))
	}
}

#Function for producing input string for RAxML and calling up the program (requires that RAxML be in your
#designated path). It reads in all .phylip files in the designated path, and runs RAxML for each in turn.
#Outgroups and mutation models can be specified either as a single string to be used for all loci or as a
#vector which needs to match the order of the reading of the phylip files (i.e., alphabetic/numeric).
RunRaxml <- function(raxmlPath=paste(getwd(),"/",sep=""),raxmlVersion="raxmlHPC",
	inputPath=paste(getwd(),"/",sep=""),mutationModel,outgroup,iterations,
	seed=sample(1:10000000,1),outputSeeds=FALSE,discard=FALSE){
	phylipFilesList <- list.files(inputPath,pattern="*.phylip",full.names=FALSE)
	seed.vec <- array()
	for(numb in 1:length(phylipFilesList)){
		inputFile=phylipFilesList[numb]
		seed=sample(1:10000000,1)
		currentSeed <- seed
		seed.vec <- c(seed.vec,currentSeed)
		outputFile <- paste(inputFile,".tre",sep="")
		if(mutationModel=="file"){
			mutationModelFile <- (read.table(paste(inputPath,"mutationModel.txt",sep="")))
			thisMutationModel <- as.character(mutationModelFile[[1]][numb])
		} else{
			thisMutationModel <- mutationModel
		}
		if(outgroup=="file"){
			outgroupFile <- (read.table(paste(inputPath,"outgroup.txt",sep="")))
			thisOutgroup <- as.character(outgroupFile[[1]][numb])
		} else {
			thisOutgroup <-outgroup
		}

 		systemCall2 <- system(paste(raxmlPath,raxmlVersion," -w ",inputPath," -s ",inputPath,inputFile," -n ",
 			outputFile," -m ",thisMutationModel," -o ",thisOutgroup," -f d -N ",iterations," -p ",
 			currentSeed,sep=""))

		#Dicard inessential RAxML output
		if(discard==TRUE){
			if(file.exists(paste(inputPath,inputFile,".reduced",sep=""))){
				unlink(paste(inputPath,inputFile,".reduced",sep=""))
			}
			if(file.exists(paste(inputPath,"RAxML_info.",inputFile,".tre",sep=""))){
				unlink(paste(inputPath,"RAxML_info.",inputFile,".tre",sep=""))
			}
			for(run in 0:(iterations-1)){
				if(file.exists(paste(inputPath,"RAxML_log.",inputFile,".tre.RUN.",run,sep=""))){
					unlink(paste(inputPath,"RAxML_log.",inputFile,".tre.RUN.",run,sep=""))
				}
				if(file.exists(paste(inputPath,"RAxML_parsimonyTree.",inputFile,".tre.RUN.",run,sep=""))){
					unlink(paste(inputPath,"RAxML_parsimonyTree.",inputFile,".tre.RUN.",run,sep=""))
				}
				if(file.exists(paste(inputPath,"RAxML_result.",inputFile,".tre.RUN.",run,sep=""))){
					unlink(paste(inputPath,"RAxML_result.",inputFile,".tre.RUN.",run,sep=""))
				}
			}
		}
	}
	if(outputSeeds==TRUE){
		write.table(seed.vec[-1],file=paste(inputPath,"RAxML.seeds.txt",sep=""))
	}
	return(systemCall2)
}

#This takes a path to .tre files and merges all trees into a single file called trees.tre
MergeTrees <- function(treesPath){	
	filenames <- list.files(treesPath,pattern="*.tre",full.names=FALSE)
	vecotr <- lapply(paste(treesPath,filenames,sep=""),read.table)
	for (treeRep in 1:length(vecotr)){
		write.table(vecotr[[treeRep]],file=paste(treesPath,"trees.tre",sep=""),append=TRUE,
			row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
}

#If only using a subset of models, this partitions migrationArrayMap to match the chosen model range
GenerateMigrationArrayMapTrunc<-function(migrationArrayMap,modelRange){
	migrationArrayMapTrunc=cbind(1:length(modelRange),migrationArrayMap[modelRange,-1])
	colnames(migrationArrayMapTrunc)<-colnames(migrationArrayMap)
	return(migrationArrayMapTrunc)
}

#This function calculates a scaled number of subsample replicates that
#is corrected based on the average original sample sizes from which subsamples
#are taken. The larger the original sample sizes, the more independent
#that subsamples are, and the more that each subsample should contribute
#to the information value in a dataset (and thus to the lnL). Inputs are popAssignments (the number
#of individuals subsampled per population), totalPopVector (the number of samples per
#locus from which subsamples are taken), and subsamplesPerGene. 
#The minumum possible scaled nreps = subsamplesPerGene,
#which results when the total sample size is equal to popAssignments[[1]] (in this case, matches among
#subsamples are simply averaged).
GetScaledNreps<-function(popAssignments,subsamplesPerGene=1,totalPopVector){
	#Get effective number of independent samples possible
	effectiveNumSamples<- (totalPopVector / length(popAssignments[[1]])) / mean(popAssignments[[1]])

	#Get effective sample size scaler (e.g., the proportion of a subsample replicate
	#that is possibly 'independent')
	eSamplesPerSubsample<-effectiveNumSamples / subsamplesPerGene

	#If the effective number of independent samples is greater than the number
	#of subsample replicates, set scalar to one (if this is the case, one
	#should probably increase the number of subsamples to allow all the
	#potentially unique subsamples to be sampled
	if(eSamplesPerSubsample > 1){
		eSamplesPerSubsample<- 1
	}

	#Calculate scaled nreps by taking the reciprocal of the effective 
	#sample scalar. This is the number by which you divide matches summed
	#across subsamples for each locus
	eNreps<- 1 / eSamplesPerSubsample

	#Return the scaled nreps (default is based on the smallest population size)
	return(eNreps)
}

#SummaryFn function (When we want to sum lnLs across subsamples and divide by the scaled number
#of replicates)
SumDivScaledNreps<-function(localVector,popAssignments,subsamplesPerGene=1,totalPopVector){
	eNreps<-GetScaledNreps(popAssignments=popAssignments,subsamplesPerGene=subsamplesPerGene,
		totalPopVector=totalPopVector)
	summaryFn<-sum(localVector) / eNreps 
	return(summaryFn)
}

#This function takes output from an exhaustive search and assembles AICs and Likelihoods for a given set
#of models into a table
ExtractAICs<-function(result=result,migrationArray=migrationArray,modelRange=c(1:length(migrationArray))){

	#Pull out the overall AICs
	AIC<-grep("objective",result,value=TRUE)
	AIC<-gsub(".+objective = ","",AIC)
	AIC<-gsub(", solution.+","",AIC)
	
	#Construct dataframe consisting of AICs and model descriptions for each model
  	overall.results<-data.frame
	params.K<-rep(NA, length(AIC))
	params.vector<-rep(NA, length(AIC))
	for (i in sequence(length(AIC))) {
		AIC[i]
		KAll(migrationArray[[i]])
		params.K[i]<-KAll(migrationArray[[i]])
		params.vector[i]<-paste(MsIndividualParameters(migrationArray[[i]]), sep=" ", collapse=" ")
	}
	models<-as.character(modelRange)
	params.K<-as.character(params.K)
	
	#Back-transform AICs to get Likelihoods
	options(digits=12)
	lnL<- as.character(round((as.numeric(AIC) / -2) + as.numeric(params.K),12))
	
	#Combine data
	overall.results<-data.frame(models, AIC, lnL, params.K, params.vector,stringsAsFactors=FALSE)
	overall.results[,1]<-as.numeric(overall.results[,1])
	overall.results[,2]<-as.numeric(overall.results[,2])
	overall.results[,3]<-as.numeric(overall.results[,3])
	overall.results[,4]<-as.numeric(overall.results[,4])

	return(overall.results)
}

#This function takes output from an exhaustive search and assembles AICs and Likelihoods for a given set
#of models into a table
ExtractGridAICs<-function(result=result,migrationArray=migrationArray,modelRange=c(1:length(migrationArray))){

	#Pull out best AIC for each model 
	AIC<-c()
	for(rep in 1:length(result)){
		AIC<-append(AIC,result[[rep]]$AIC[1])
	}
	
	#Construct dataframe consisting of AICs and model descriptions for each model
  	overall.results<-data.frame
	params.K<-rep(NA, length(AIC))
	params.vector<-rep(NA, length(AIC))
	for (i in sequence(length(AIC))) {
		AIC[i]
		KAll(migrationArray[[i]])
		params.K[i]<-KAll(migrationArray[[i]])
		params.vector[i]<-paste(MsIndividualParameters(migrationArray[[i]]), sep=" ", collapse=" ")
		params.vector[i]<-gsub("n0multiplier_1","",params.vector[i]) #toss first n0multiplier parameter
	}
	models<-as.character(modelRange)
	params.K<-as.character(params.K)
	
	#Back-transform AICs to get Likelihoods
	options(digits=12)
	lnL<- as.character(round((as.numeric(AIC) / -2) + as.numeric(params.K),12))
	
	#Combine data
	overall.results<-data.frame(models, AIC, lnL, params.K, params.vector,stringsAsFactors=FALSE)
	overall.results[,1]<-as.numeric(overall.results[,1])
	overall.results[,2]<-as.numeric(overall.results[,2])
	overall.results[,3]<-as.numeric(overall.results[,3])
	overall.results[,4]<-as.numeric(overall.results[,4])

	return(overall.results)
}

#This function takes output from an exhaustive search and assembles parameter indexes and estimates
#based on a set of models. A list containing two tables is outputted: one containing parameter indexes
#and one containing parameter estimates
#This function takes output from an exhaustive search and assembles parameter indexes and estimates
#based on a set of models. A list containing two tables is outputted: one containing parameter indexes
#and one containing parameter estimates
ExtractParameters<-function(migrationArray=migrationArray,result=result,popVector){

	############MAKE COLUMN HEADERS FOR MIGRATION BASED ON FULL SQUARE MATRICES 
	############(INCLUDING ALL BUT THE DIAGONAL NA's)
	
	npop<-length(popVector)
	#Establish column names for these parameters (keep square matrices method)
	allMigs<-list()
	for(thisMatrix in 1:(npop - 1)){ #for each possible matrix, given npop
		currentMigs<-array(NA,dim=c(npop,npop)) # build dimensions of the current migration
		allMigs[[length(allMigs)+1]] <- currentMigs #add it to the master list
	}
	#Fill in an array with migration parameter headers
	migColnames=c()
	for(thisMatrix in 1:length(allMigs)){ #for each migration matrix in the model
		presentMatrix<-allMigs[[1]]
		for(migFrom in 1:length(presentMatrix[,1])){ #for each row
			for(migTo in 1:length(presentMatrix[1,])){ #for each column
				if(migFrom !=migTo){ #skip diagonal migrations
					migColnames<-append(migColnames,paste("m",thisMatrix,"_",migFrom,".",migTo,sep=""))
				}
			}
		}
	}	
	
	
	############MAKE COLUMN HEADERS FOR COLLAPSE AND N0MULTI BASED ON FULLY RESOLVED TREE
	#First, make rectangular matrix to match fully resolved case 
	allCollapses<-array(NA,dim=c(npop,npop - 1))
	
	#Fill in arrays with parameter headers
	collapseColnames=c()
	n0multiColnames=c()
	for(cols in 1:length(allCollapses[1,])){ #for each column
		for(rows in 1:length(allCollapses[,1])){ #for each row
			collapseColnames<-append(collapseColnames,paste("tau_t",cols,".",rows,sep=""))
			n0multiColnames<-append(n0multiColnames,paste("n_t",cols,".",rows,sep=""))
		}
	}	
	
	
	#############FILL IN A MATRIX WITH THE MIGRATION PARAMETER INDEXES
	#First establish the total number of columns required to fit the migration parameters
	migParmsNcol<-(npop^2 - npop) * (npop - 1)
	
	#Create empty matrix to fill
	migParms<-data.frame(matrix(NA,nrow=length(migrationArray),ncol=migParmsNcol)) #for parameter indexes
	currentModelCount=0
	
	for(currentModel in 1:length(migrationArray)){ #loop through each model
		currentModelCount=currentModelCount+1
		counter<-0
		for(thisMatrix in 1:length(migrationArray[[currentModel]]$migrationArray[1,1,])){ #loop through each 
			#migration matrix in the model. Note that the first two numbers give rows and columns of a 
			#matrix, the last number gives the number of matrices (current and historical)
		currentMatrix<-migrationArray[[currentModel]]$migrationArray[,,thisMatrix] #rename current matrix
			for(migFrom in 1:length(currentMatrix[,1])){ #for each row	
				for(migTo in 1:length(currentMatrix[1,])){ #for each column		
					if(migFrom != migTo){ #skip migrations into same population	
						counter<-counter+1 #keep track of column to deposit values	
						migParms[currentModelCount,counter]<-currentMatrix[migFrom,migTo] #pull out the migration indexes
					}
				}
			}
		}
	}			
	colnames(migParms)<-migColnames
			

	#############FILL IN A MATRIX WITH THE COLLAPSE AND N0MULTI PARAMETER INDEXES
	
	#Create empty matrices to fill
	collapseParms<-data.frame(matrix(NA,nrow=length(migrationArray),ncol=(npop * (npop - 1)))) 
	n0multiParms<-data.frame(matrix(NA,nrow=length(migrationArray),ncol=(npop * (npop - 1))))
	currentModelCount=0
	
	for(currentModel in 1:length(migrationArray)){ #loop through each model
		collapseCurrentMatrix<-migrationArray[[currentModel]]$collapse #current model collpases
		n0multiCurrentMatrix<-migrationArray[[currentModel]]$n0multiplier #current model n0multis
		currentModelCount=currentModelCount+1
		counter<-0
		for(cols in 1:length(collapseCurrentMatrix[1,])){ #for each column
			for(rows in 1:length(collapseCurrentMatrix[,1])){ #for each row	
				counter<-counter+1 #keep track of column to deposit values	
				collapseParms[currentModelCount,counter]<-collapseCurrentMatrix[rows,cols] #pull out indexes
				n0multiParms[currentModelCount,counter]<-n0multiCurrentMatrix[rows,cols] #pull out indexes
			}
		}
	}			
	colnames(collapseParms)<-collapseColnames
	colnames(n0multiParms)<-n0multiColnames
	

	
	##########################MAKE MATRICES FILLED WITH THE PARAMETER VALUES, BASED ON THE INDEXES
	#Create data.frames to substitute in with the parameter values
	collapseParmsValues<-collapseParms
	n0multiParmsValues<-n0multiParms
	migParmsValues<-migParms
	
	#Now loop through each  (i.e., each row in Parms matrices)
	for(eachCol in 1:nrow(collapseParms)){ 
	
		#First, make matrix of the collapse estimates
		currentSol<-exp(result[[eachCol]]$solution) #back-transform the logged parameter estimates
		if(length(currentSol)> 0){ #if there are parameter estimates left
				
			#Because different collapse parameters are denoted by different columns rather than by different index values,
			#for each new time period (i.e., column), if there is a collapse, create a new index number so that a distinct
			#tau is assigned to each coalescent event.
			makeIndexes<-as.numeric(collapseParms[eachCol,])
			count1<-1
			count2<-npop
			count3<-1
			for(rep in 1:(length(makeIndexes) / npop)){
				if(length(makeIndexes[count1:count2][which(makeIndexes[count1:count2]==1)]) > 0){
					makeIndexes[count1:count2][which(makeIndexes[count1:count2]==1)] <- count3
					count3<-count3 + 1
				}
				count1<-count1 + npop
				count2<-count2 + npop
			}

			#Now calulate the number of unique parameters
			uniqueCollapse<-makeIndexes[!is.na(makeIndexes)]
			uniqueCollapse<-uniqueCollapse[which(uniqueCollapse > 0)]
			uniqueCollapse<-length(unique(uniqueCollapse))
			if(uniqueCollapse > 0){ #If there are any unique parameters
				for(rep1 in 1:uniqueCollapse){ #For each unique parameter
					for(rep2 in sequence(1:length(makeIndexes))){ #Go through each parameter position
						if(!is.na(makeIndexes[rep2])){
							if(makeIndexes[rep2]==rep1){ #if the index matches the unique parameter
								collapseParmsValues[eachCol,rep2]<-currentSol[rep1] #change the index to the corresponding parameter value
							}
						}
					}
				}
			}
		}

		#Now take away from the parameters vector the used up parameter values
		currentSol<-tail(currentSol,length(currentSol) - uniqueCollapse)
			
	
		#Second, make matrix of the n0multiplier estimates
		if(length(currentSol)> 0){ #if there are parameter estimates left
			#Calulate the number of unique parameters
			uniqueN0multi<-as.numeric(n0multiParms[eachCol,]) 
			uniqueN0multi<-uniqueN0multi[!is.na(uniqueN0multi)]
			uniqueN0multi<-uniqueN0multi[which(uniqueN0multi > 0)]
			uniqueN0multi<-length(unique(uniqueN0multi))
			if(uniqueN0multi > 0){ #If there are any unique parameters
				for(rep1 in 1:uniqueN0multi){ #For each unique parameter
					for(rep2 in sequence(1:ncol(n0multiParms))){ #Go through each parameter position
						if(!is.na(as.numeric(n0multiParms[eachCol,rep2]))){
							if(as.numeric(n0multiParms[eachCol,rep2])==rep1){ #if the index matches the unique parameter
								n0multiParmsValues[eachCol,rep2]<-currentSol[rep1] #change the index to the corresponding parameter value
							}
						}
					}
				}
			}
		}
		
		#Now take away from the parameters vector the used up parameter values
		currentSol<-tail(currentSol,length(currentSol) - uniqueN0multi)
	
	
		#Third, make matrix of the migration estimates
		if(length(currentSol)> 0){ #if there are parameter estimates left
			#Calulate the number of unique parameters
			uniqueMig<-as.numeric(migParms[eachCol,]) 
			uniqueMig<-uniqueMig[!is.na(uniqueMig)]
			uniqueMig<-uniqueMig[which(uniqueMig > 0)]
			uniqueMig<-length(unique(uniqueMig))
			if(uniqueMig > 0){ #If there are any unique parameters
				for(rep1 in 1:uniqueMig){ #For each unique parameter
					for(rep2 in sequence(1:ncol(migParms))){ #Go through each parameter position
						if(!is.na(as.numeric(migParms[eachCol,rep2]))){
							if(as.numeric(migParms[eachCol,rep2])==rep1){ #if the index matches the unique parameter
								migParmsValues[eachCol,rep2]<-currentSol[rep1] #change the index to the corresponding parameter value
							}
						}
					}
				}
			}
		}
	}
	
	
	##########################CBIND ALL PARAMETER INDEX COLUMNS TOGETHER
	parameterIndexes<-cbind(collapseParms,n0multiParms,migParms)
	#Add an "I" to colnames for the indexes (to distinguish them from the parameter values colnames)
	for(rep in 1:ncol(parameterIndexes)){
		colnames(parameterIndexes)[rep]<-paste(colnames(parameterIndexes)[rep],"_I",sep="")
	}
	parameterValues<-cbind(collapseParmsValues,n0multiParmsValues,migParmsValues)
	parameterAll<-cbind(parameterValues,parameterIndexes)
	
	return(list(parameterValues,parameterIndexes))
}



#This function takes output from an exhaustive search and assembles parameter indexes and estimates
#based on a set of models. A list containing two tables is outputted: one containing parameter indexes
#and one containing parameter estimates
ExtractGridParameters<-function(migrationArray=migrationArray,result=result,popVector){

	############MAKE COLUMN HEADERS FOR MIGRATION BASED ON FULL SQUARE MATRICES 
	############(INCLUDING ALL BUT THE DIAGONAL NA's)
	
	npop<-length(popVector)
	#Establish column names for these parameters (keep square matrices method)
	allMigs<-list()
	for(thisMatrix in 1:(npop - 1)){ #for each possible matrix, given npop
		currentMigs<-array(NA,dim=c(npop,npop)) # build dimensions of the current migration
		allMigs[[length(allMigs)+1]] <- currentMigs #add it to the master list
	}
	#Fill in an array with migration parameter headers
	migColnames=c()
	for(thisMatrix in 1:length(allMigs)){ #for each migration matrix in the model
		presentMatrix<-allMigs[[1]]
		for(migFrom in 1:length(presentMatrix[,1])){ #for each row
			for(migTo in 1:length(presentMatrix[1,])){ #for each column
				if(migFrom !=migTo){ #skip diagonal migrations
					migColnames<-append(migColnames,paste("m",thisMatrix,"_",migFrom,".",migTo,sep=""))
				}
			}
		}
	}	
	
	
	############MAKE COLUMN HEADERS FOR COLLAPSE BASED ON FULLY RESOLVED TREE
	#First, make rectangular matrix to match fully resolved case 
	allCollapses<-array(NA,dim=c(npop,npop - 1))
	
	#Fill in arrays with parameter headers
	collapseColnames=c()
	for(cols in 1:length(allCollapses[1,])){ #for each column
		for(rows in 1:length(allCollapses[,1])){ #for each row
			collapseColnames<-append(collapseColnames,paste("tau_t",cols,".",rows,sep=""))
		}
	}	
	
	
	#############FILL IN A MATRIX WITH THE MIGRATION PARAMETER INDEXES
	#First establish the total number of columns required to fit the migration parameters
	migParmsNcol<-(npop^2 - npop) * (npop - 1)
	
	#Create empty matrix to fill
	migParms<-data.frame(matrix(NA,nrow=length(migrationArray),ncol=migParmsNcol)) #for parameter indexes
	currentModelCount=0
	
	for(currentModel in 1:length(migrationArray)){ #loop through each model
		currentModelCount=currentModelCount+1
		counter<-0
		for(thisMatrix in 1:length(migrationArray[[currentModel]]$migrationArray[1,1,])){ #loop through each 
			#migration matrix in the model. Note that the first two numbers give rows and columns of a 
			#matrix, the last number gives the number of matrices (current and historical)
		currentMatrix<-migrationArray[[currentModel]]$migrationArray[,,thisMatrix] #rename current matrix
			for(migFrom in 1:length(currentMatrix[,1])){ #for each row	
				for(migTo in 1:length(currentMatrix[1,])){ #for each column		
					if(migFrom != migTo){ #skip migrations into same population	
						counter<-counter+1 #keep track of column to deposit values	
						migParms[currentModelCount,counter]<-currentMatrix[migFrom,migTo] #pull out the migration indexes
					}
				}
			}
		}
	}			
	colnames(migParms)<-migColnames
			

	#############FILL IN A MATRIX WITH THE COLLAPSE PARAMETER INDEXES
	
	#Create empty matrices to fill
	collapseParms<-data.frame(matrix(NA,nrow=length(migrationArray),ncol=(npop * (npop - 1)))) 
	currentModelCount=0
	
	for(currentModel in 1:length(migrationArray)){ #loop through each model
		collapseCurrentMatrix<-migrationArray[[currentModel]]$collapse #current model collpases
		currentModelCount=currentModelCount+1
		counter<-0
		for(cols in 1:length(collapseCurrentMatrix[1,])){ #for each column
			for(rows in 1:length(collapseCurrentMatrix[,1])){ #for each row	
				counter<-counter+1 #keep track of column to deposit values	
				collapseParms[currentModelCount,counter]<-collapseCurrentMatrix[rows,cols] #pull out indexes
			}
		}
	}			
	colnames(collapseParms)<-collapseColnames
	

	
	##########################MAKE MATRICES FILLED WITH THE PARAMETER VALUES, BASED ON THE INDEXES
	#Create data.frames to substitute in with the parameter values
	collapseParmsValues<-collapseParms
	migParmsValues<-migParms
	
	#Now loop through each  (i.e., each row in Parms matrices)
	for(eachCol in 1:nrow(collapseParms)){ 

		#First, make matrix of the collapse estimates
		#Make vector of top parameter estimates from the grid (based on the mean of the top five estimates for each parm)
		currentSol<-c()
		positionOfNonParameters <- grep("lnL", colnames(result[[eachCol]]))[1] #find slopes
		if(!is.na(positionOfNonParameters)){
			currentResult<-result[[eachCol]][2:(positionOfNonParameters - 1)] #toss slopes
		}else{
			currentResult<-result[[eachCol]][2:length(result[[eachCol]])] #if only 1 popAssignment used
		}

		for(rep in 1:ncol(currentResult)){
			currentSol<-append(currentSol,mean(currentResult[1:5,rep])) #Note that grid prams are already back-transformed (i.e., logged)
		}
		if(length(currentSol)> 0){ #if there are parameter estimates left
				
			#Because different collapse parameters are denoted by different columns rather than by different index values,
			#for each new time period (i.e., column), if there is a collapse, create a new index number so that a distinct
			#tau is assigned to each coalescent event.
			makeIndexes<-as.numeric(collapseParms[eachCol,])
			count1<-1
			count2<-npop
			count3<-1
			for(rep in 1:(length(makeIndexes) / npop)){
				if(length(makeIndexes[count1:count2][which(makeIndexes[count1:count2]==1)]) > 0){
					makeIndexes[count1:count2][which(makeIndexes[count1:count2]==1)] <- count3
					count3<-count3 + 1
				}
				count1<-count1 + npop
				count2<-count2 + npop
			}

			#Now calulate the number of unique parameters
			uniqueCollapse<-makeIndexes[!is.na(makeIndexes)]
			uniqueCollapse<-uniqueCollapse[which(uniqueCollapse > 0)]
			uniqueCollapse<-length(unique(uniqueCollapse))
			if(uniqueCollapse > 0){ #If there are any unique parameters
				for(rep1 in 1:uniqueCollapse){ #For each unique parameter
					for(rep2 in sequence(1:length(makeIndexes))){ #Go through each parameter position
						if(!is.na(makeIndexes[rep2])){
							if(makeIndexes[rep2]==rep1){ #if the index matches the unique parameter
								collapseParmsValues[eachCol,rep2]<-currentSol[rep1] #change the index to the corresponding parameter value
							}
						}
					}
				}
			}
		}

		#Now take away from the parameters vector the used up parameter values
		currentSol<-tail(currentSol,length(currentSol) - uniqueCollapse)	
	
		#Second, make matrix of the migration estimates
		if(length(currentSol)> 0){ #if there are parameter estimates left
			#Calulate the number of unique parameters
			uniqueMig<-as.numeric(migParms[eachCol,]) 
			uniqueMig<-uniqueMig[!is.na(uniqueMig)]
			uniqueMig<-uniqueMig[which(uniqueMig > 0)]
			uniqueMig<-length(unique(uniqueMig))
			if(uniqueMig > 0){ #If there are any unique parameters
				for(rep1 in 1:uniqueMig){ #For each unique parameter
					for(rep2 in sequence(1:ncol(migParms))){ #Go through each parameter position
						if(!is.na(as.numeric(migParms[eachCol,rep2]))){
							if(as.numeric(migParms[eachCol,rep2])==rep1){ #if the index matches the unique parameter
								migParmsValues[eachCol,rep2]<-currentSol[rep1] #change the index to the corresponding parameter value
							}
						}
					}
				}
			}
		}
	}
	
	
	##########################CBIND ALL PARAMETER INDEX COLUMNS TOGETHER
	parameterIndexes<-cbind(collapseParms,migParms)
	#Add an "I" to colnames for the indexes (to distinguish them from the parameter values colnames)
	for(rep in 1:ncol(parameterIndexes)){
		colnames(parameterIndexes)[rep]<-paste(colnames(parameterIndexes)[rep],"_I",sep="")
	}
	parameterValues<-cbind(collapseParmsValues,migParmsValues)
	parameterAll<-cbind(parameterValues,parameterIndexes)
	
	return(list(parameterValues,parameterIndexes))
}





###############################Post-analysis Functions####################

#This concatenates phrapl results, and also adds to these dAIC, wAIC, and model ranks for a set of models
ConcatenateResults<-function(rdaFilesPath="./",rdaFiles=NULL,outFile=NULL,addAICweights=TRUE,
		rmNaParameters=TRUE,addTime.elapsed=FALSE){
	#If a vector of rda file names is not provided, then read them in from the given path
	if(is.null(rdaFiles)){
		rdaFiles<-grep(".rda",list.files(rdaFilesPath),value=TRUE)
	}
	totalData<-data.frame()
	result <- c()
	time.elapsed <- c()
	for (rep in 1:length(rdaFiles)){
		load(paste(rdaFilesPath,rdaFiles[rep],sep="")) #Read model objects

		#Combine results from the rda into dataframe
		overall.results<-result[[2]][order(result[[2]][,1]),] #sort by model number (same order as parameters)
		if(addTime.elapsed==TRUE){
			current.results<-cbind(overall.results,elapsedHrs=time.elapsed[,2],result[[3]],result[[4]])
		}else{
			current.results<-cbind(overall.results,result[[3]],result[[4]])
		}
		totalData<-rbind(totalData,current.results) #Add current.results to totalData
		totalData<-totalData[order(totalData$models),] #Make sure totalData is all sorted by model
	}		
	#Add model ranks, dAIC, and wAIC
	if(addAICweights==TRUE){
		dAICRankWeights<-GetAICweights(totalData)
		totalData<-cbind(totalData[,1:3],dAICRankWeights,totalData[,4:ncol(totalData)])					
	}
		
	#Remove Parameters that only have NAs for values
	if(rmNaParameters==TRUE){
		totalData<-RemoveParameterNAs(totalData)
	}
	
	#Print table
	if(!is.null(outFile)){
		write.table(totalData[order(totalData$dAIC,totalData$models),],
			outFile,quote=FALSE,sep="\t",row.names=FALSE)
	}
	
	return(totalData[order(totalData$dAIC,totalData$models),])
}
	            
	
#Function used for calculating AIC weights
getExponent<-function(x){
	return(exp(-0.5 * x))
	}

#Remove parameter columns that contain only NAs
RemoveParameterNAs<-function(totalData){
	columnsLeft=ncol(totalData) #update number of remaining columns
	continue=TRUE #Not to last column?
	proceed=TRUE
	colnum=0 #Present column
	while(continue){
		if(proceed==TRUE){
			colnum<-colnum + 1
		}
		NAlist<-is.na(totalData[,colnum])
		proceed<-TRUE
		if(length(NAlist[which(NAlist==TRUE)])==length(NAlist)){ #if there are only NAs in column
			
			totalData<-totalData[,-colnum] #prune column from dataset
			columnsLeft<-ncol(totalData)
			colnum<-colnum #don't proceed to increase current column (because just tossed one)
			proceed=FALSE
	
		}

		if(colnum >= columnsLeft){
			continue=FALSE
		}
		
		if(continue==FALSE){
			NAlist<-is.na(totalData[,length(totalData)])
			if(length(NAlist[which(NAlist==TRUE)])==length(NAlist)){
				totalData<-totalData[,-length(totalData)]
			}
		}
	}
	
	return(totalData)
}

#Calculate dAICs, model ranks, and model weights
GetAICweights<-function(totalData=totalData){
	#Calculate change in AIC values from best model
	totalData<-totalData[order(totalData$AIC,totalData$models),] #sort by AIC, then by model
	dAIC=array(0)
	for(rep4 in 2:nrow(totalData)){
		dAIC<-c(dAIC,round((totalData$AIC[rep4] - totalData$AIC[1]),digits=3))
	}
	
	#Add a model rank column				
	rank=array(1)
	for(rep3 in 2:nrow(totalData)){
		if(totalData$AIC[rep3]==totalData$AIC[rep3-1]){
			rank=c(rank,rank[rep3-1])
		}else{
			rank=c(rank,rank[rep3-1]+1)
		}
	}

	#Calculate Weights across subsamples
	currentExp<-sapply(dAIC,getExponent)
	wAIC<-currentExp / sum(currentExp)
	dAICRankWeights<-data.frame(cbind(models=totalData$models,rank,dAIC,wAIC))
	dAICRankWeights<-dAICRankWeights[order(dAICRankWeights$models),]
	return(dAICRankWeights[,-1])
}

#Get weights for three model types (Isolation only, migration only, isolation + migration)
#based on AIC weights. Output can either be "indexes" or "weights" (where the indexes are logical
#columns for the three model types).
GetTrianglePlotWeights<-function(totalData=totalData,output="weights"){	
	#Get datasets containing isolation AND migration parameters
	currentIsoMig<-grepl(".*collapse.*migration.*",totalData$params.vector)

	#Get datasets containing ONLY isolation or ONLY migration
	currentIsoAll<-grepl(".*collapse.*",totalData$params.vector)
	currentMigAll<-grepl(".*migration.*",totalData$params.vector)
	currentIso<-array()
	currentMig<-array()
	for(rep in 1:length(totalData$params.vector)){
		if(currentIsoAll[rep]==TRUE && currentMigAll[rep]==FALSE){
			currentIso[rep]<-TRUE
		}else{
			currentIso[rep]<-FALSE
		}
		
		if(currentMigAll[rep]==TRUE && currentIsoAll[rep]==FALSE){
			currentMig[rep]<-TRUE
		}else{
			currentMig[rep]<-FALSE
		}
	}
	modelIndexes<-data.frame(cbind(currentIso,currentMig,currentIsoMig))
	colnames(modelIndexes)<-c("iso","mig","isoANDmig")
	
	#calculate model type weights by summing weights under the three types of weighted models
	if(output=="weights"){
		weightIso<-sum(totalData$wAIC[which(modelIndexes$iso==TRUE)])
		weightMig<-sum(totalData$wAIC[which(modelIndexes$mig==TRUE)])
		weightIsoMig<-sum(totalData$wAIC[which(modelIndexes$isoANDmig==TRUE)])
		plotWeights<-data.frame(cbind(weightIso,weightMig,weightIsoMig))
		return(plotWeights) #same as old currentPlotWeights
	}else{
		return(modelIndexes)
	}	
}

#Calculate model averaged parameter values across a set of models
CalculateModelAverages<-function(totalData=totalData,parmStartCol=10){
	#Add model averaged parameter estimates to a dataframe
	totalData<-totalData[order(totalData$AIC,totalData$models),] #sort parameters by AIC to match totalData
	totalData<-totalData[,grep(".*_I",colnames(totalData),invert=TRUE)] #Remove parameter index columns
										
	#Remove parameter columns that contain only NAs
	totalData<-RemoveParameterNAs(totalData)
		
	#Calculate model averaged parameter estimates (NA's are ignored if present in some models)					
	modelAverages<-na.pass(totalData[,parmStartCol:ncol(totalData)] * totalData$wAIC)
	modelAverages<-apply(modelAverages,2,sum,na.rm=TRUE)
					
	#Convert modelAverages into a dataframe
	modelAveragesMat<-data.frame(matrix(modelAverages,nrow=1,ncol=length(names(modelAverages))))
	colnames(modelAveragesMat)<-names(modelAverages)
	modelAverages<-modelAveragesMat
	return(modelAverages)
}


################Functions summarizing across subsamples, models, and treatments

#This function summarizes (mean, median, and standard deviation) triangle plot weights and model-averaged parameter
#estimates across subsamples of the same treatment. treatmentVec gives a vector of treatments assigned to each model. 
#paramColVec gives the column numbers that are to be summarized with mean, median, and standard deviation.
SummarizeAveragedModelsAcrossSubsamples<-function(plotWeights=plotWeights,paramColVec=c(7:ncol(plotWeights)),treatmentVec=NULL){
	
	#Define all unique treatments (if any)
	if(!is.null(treatmentVec)){
		utreat=unique(treatmentVec)
	}else{
		treatmentVec<-rep("treat",nrow(plotWeights))
		utreat<-unique(treatmentVec)
	}

	#Calculate mean, median, and sd across subsamples for each treatment
	medianWeights<-data.frame()
	meanWeights<-data.frame()
	sdWeights<-data.frame()

	for(rep10 in 1:length(utreat)){
		currentPlotWeights<-plotWeights[which(treatmentVec==utreat[rep10]),paramColVec]
		currentDescriptors<-plotWeights[which(treatmentVec==utreat[rep10]),1:(paramColVec[1] - 1)]
	
		currentMedianWeights<-sapply(currentPlotWeights,median,na.rm=TRUE)
		currentMedianWeightsDF<-data.frame(matrix(currentMedianWeights,nrow=1,ncol=length(names(currentMedianWeights))))
		colnames(currentMedianWeightsDF)<-names(currentMedianWeights)
		currentMedianWeights<-cbind(currentDescriptors[1,],currentMedianWeightsDF)
		medianWeights<-rbind(medianWeights,currentMedianWeights)
	
		currentMeanWeights<-sapply(currentPlotWeights,mean,na.rm=TRUE)
		currentMeanWeightsDF<-data.frame(matrix(currentMeanWeights,nrow=1,ncol=length(names(currentMeanWeights))))
		colnames(currentMeanWeightsDF)<-names(currentMeanWeights)
		currentMeanWeights<-cbind(currentDescriptors[1,],currentMeanWeightsDF)
		meanWeights<-rbind(meanWeights,currentMeanWeights)

		currentSdWeights<-sapply(currentPlotWeights,sd,na.rm=TRUE)
		currentSdWeightsDF<-data.frame(matrix(currentSdWeights,nrow=1,ncol=length(names(currentSdWeights))))
		colnames(currentSdWeightsDF)<-names(currentSdWeights)
		currentSdWeights<-cbind(currentDescriptors[1,],currentSdWeightsDF)
		sdWeights<-rbind(sdWeights,currentSdWeights)
	}
	summaryWeights<-list(medianWeights,meanWeights,sdWeights)
}

#This function summarizes (e.g., mean, median, and standard deviation) parameter values for each model across subsamples
#treatmentVec gives a vector of treatments assigned to each model. paramColVec gives the column numbers that are to be
#summarized with mean, median, and standard deviation. timeColVec gives the column numbers that are to be summarized by 
#taking the sum across subsamples (e.g., for elapsed time per model). descriptionColVec gives the column numbers for
#columns that should be included in the new table, but which aren't to be summarized. This function also adusts the model 
#rank column based on the mean and median AIC weight across subsamples. Output for this function is a list of three
#dataframes for mean, median, and standard deviation values.
SummarizeParametersAcrossSubsamples<-function(totalData,treatmentVec=NULL,paramColVec=NULL,
	timeColVec=c(15),descriptionColVec=c(1:7,14)){
	if(is.null(paramColVec)) {
		paramColVec<-c(8:13,16:ncol(totalData))
	}
	#Define all unique treatments (if any)
	if(!is.null(treatmentVec)){
		utreat=unique(treatmentVec)
	}else{
		treatmentVec<-rep("treatment",nrow(totalData))
		utreat<-unique(treatmentVec)
	}

	#Calculate mean, median, and sd across subsamples for each model within each treatment
	medianModel<-data.frame()
	meanModel<-data.frame()
	sdModel<-data.frame()

	for(rep11 in 1:length(utreat)){ #for each treatment
		currentData<-totalData[which(treatmentVec==utreat[rep11]),]
		for(rep12 in 1:length(unique(currentData$models))){ #for each model
			currentModel<-currentData[which(currentData$models==unique(currentData$models)[rep12]),]
		
			#isolate columns to average
			currentModelAv<-currentModel[,paramColVec]		

			#sum up elapsed time for each model
			elapsedTime<-sapply(currentModel[,timeColVec],sum,na.rm=TRUE)
			elapsedTimeDF<-data.frame(matrix(elapsedTime,nrow=1,ncol=length(names(elapsedTime))))
			colnames(elapsedTimeDF)<-names(elapsedTime)		
		
			#Get median
			currentMedianModel<-sapply(currentModelAv,median,na.rm=TRUE)
			currentMedianModelDF<-data.frame(matrix(currentMedianModel,nrow=1,ncol=length(names(currentMedianModel))))
			colnames(currentMedianModelDF)<-names(currentMedianModel)
			currentMedianModel<-cbind(currentModel[1,descriptionColVec],elapsedTimeDF,currentMedianModelDF)
			medianModel<-rbind(medianModel,currentMedianModel)

			#Get mean
			currentMeanModel<-sapply(currentModelAv,mean,na.rm=TRUE)
			currentMeanModelDF<-data.frame(matrix(currentMeanModel,nrow=1,ncol=length(names(currentMeanModel))))
			colnames(currentMeanModelDF)<-names(currentMeanModel)
			currentMeanModel<-cbind(currentModel[1,descriptionColVec],elapsedTimeDF,currentMeanModelDF)
			meanModel<-rbind(meanModel,currentMeanModel)
		
			#Get sd
			currentSdModel<-sapply(currentModelAv,sd,na.rm=TRUE)
			currentSdModelDF<-data.frame(matrix(currentSdModel,nrow=1,ncol=length(names(currentSdModel))))
			colnames(currentSdModelDF)<-names(currentSdModel)
			currentSdModel<-cbind(currentModel[1,descriptionColVec],elapsedTimeDF,currentSdModelDF)
			sdModel<-rbind(sdModel,currentSdModel)
		}
	}
	
	#########Add rank columns to median and mean weights (based on model Averaged AIC weights)########
	#If treatmentVec=NULL, re-define utreat
	if(length(utreat)==1){
		treatmentVec<-rep("treatment",nrow(meanModel))
	}

	#Fill in ranks
	rankColNum<-match("rank",names(meanModel)) #get the column number where rank is stored
	meanModelRankSort<-data.frame()
	medianModelRankSort<-data.frame()
	for(rep13 in 1:length(utreat)){
		currentMeanModelRank<-meanModel[which(treatmentVec==utreat[rep13]),]
		currentMeanModelRankSort<-currentMeanModelRank[rev(order(currentMeanModelRank$wAIC)),]
	
		currentMedianModelRank<-medianModel[which(treatmentVec==utreat[rep13]),]
		currentMedianModelRankSort<-currentMedianModelRank[rev(order(currentMedianModelRank$wAIC)),]

		#Calculate ranks based on averaged AIC weights
		#Mean
		rank=array(1)
		for(rep14 in 2:nrow(currentMeanModelRankSort)){
			if(currentMeanModelRankSort$wAIC[rep14]==currentMeanModelRankSort$wAIC[rep14 - 1]){
				rank=rbind(rank,rank[rep14 - 1])
			}else{
				rank=rbind(rank,rank[rep14 - 1] + 1)
			}
		}	
		currentMeanModelRankSort[,rankColNum]<-rank
		meanModelRankSort<-rbind(meanModelRankSort,currentMeanModelRankSort)
	
		#Median
		rank=array(1)
		for(rep14 in 2:nrow(currentMedianModelRankSort)){
			if(currentMedianModelRankSort$wAIC[rep14]==currentMedianModelRankSort$wAIC[rep14 - 1]){
				rank=rbind(rank,rank[rep14 - 1])
			}else{
				rank=rbind(rank,rank[rep14 - 1] + 1)
			}
		}	
		currentMedianModelRankSort[,rankColNum]<-rank
		medianModelRankSort<-rbind(medianModelRankSort,currentMedianModelRankSort)
	}

	meanModel<-meanModelRankSort
	medianModel<-medianModelRankSort
	
	summaryModels<-list(medianModel,meanModel,sdModel)
}

#This function takes model AIC weights for three major model types (isolation only, migration only, 
#and isolation + migration) and adjusts these weights based on the assumption that the null expected weights
#(i.e., equal weights, scaled by the number of parameters in each model and the number of model types 
#represented in the given model set) are in the center of the triangle plot.The plot is centered around this
#null point and a correction is applied to the observed model type weights. Inputs are 1) the model weights file from phrapl, 2) the main results file from phrapl (which contains model K values), and 3) the column numbers in the model
#weights file that contain the three model type weights. Outputted are the adjusted model weights (where null
#weights are added in a new row). These weights can be used in a triangle plot
CalculateNullWeights<-function(modelWeights=NULL,modelKs=NULL,weightsColVec=c(7:9)){

##First, add datapoint based on null expectations (based only on K, and the proportion of
##each model type in the list of tested models)

	#Extract a single dataset from the total Dataframe
	#Calculate number of different models
	nModels<-length(unique(modelKs$models))
	sampleDataset<-modelKs[1:nModels,]
	
	#Make AICs based on 2 * number of parameters
	nullkAIC<-(2 * sampleDataset$params.K)
	nullkDiffAIC<-(nullkAIC - min(nullkAIC))
	nullkExpCol<-sapply(nullkDiffAIC,getExponent)
	nullkWeights<-nullkExpCol / sum(nullkExpCol)
	sampleDataset<-cbind(sampleDataset,nullkWeights)

	#calculate triangle plot weights for the current dataset by summing weights under the 
	#three types of weighted models
	modelIndexes<-GetTrianglePlotWeights(totalData=sampleDataset,output="indexes") #get model type indexes
	sampleDataset<-cbind(sampleDataset,modelIndexes) #concatenate these to sample dataset
	weightIsoNullk<-sum(sampleDataset$nullkWeights[which(sampleDataset$iso==1)])
	weightMigNullk<-sum(sampleDataset$nullkWeights[which(sampleDataset$mig==1)])
	weightIsoMigNullk<-sum(sampleDataset$nullkWeights[which(sampleDataset$isoANDmig==1)])

	#Add plot weights to a dataframe
	nullkPlotWeights<-cbind(weightIsoNullk,weightMigNullk,weightIsoMigNullk)
	nullkPlotWeightsDF<-modelWeights[NA,]
	nullkPlotWeightsDF<-nullkPlotWeightsDF[1,]
	nullkPlotWeightsDF[,weightsColVec]<-nullkPlotWeights

	#Add these null datapoints to the plotWeights dataframe
	modelWeights<-rbind(modelWeights,nullkPlotWeightsDF) 

	
	##Second, scale model weights to the centered null weights

	#Now scale the triangle plot such that the nullkPlotWeights is in the center
	nullkDev<-(1 / 3) / nullkPlotWeights[,1:3] #calculate proportional deviation of null from plot center
	modelWeightsScaled<-modelWeights[0,weightsColVec] #create empty dataframe
	for(datasets in 1:nrow(modelWeights)){ #scale old weights to proportional deviation
		modelWeightsScaled[datasets,]<-(modelWeights[datasets,weightsColVec] * nullkDev) / sum(modelWeights[datasets,
			weightsColVec] * nullkDev)
	}
	modelWeightsScaled<-cbind(modelWeights[,1:(weightsColVec[1] - 1)],modelWeightsScaled,
		modelWeights[,(weightsColVec[3] +1):ncol(modelWeights)])

return(modelWeightsScaled)
}
