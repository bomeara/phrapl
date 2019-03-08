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
	return(length(which(x[,dim(x)[2]]>=val)))
}

ColCountLastZero<-function(x,val) {
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

N0multiplierGrowthindividual<-function(collapseMatrix,complete=FALSE,n0multiplierMap,growthMap) {
	tI<-list(collapseMatrix=collapseMatrix,complete=complete,n0multiplierMap=n0multiplierMap,growthMap=growthMap)
	class(tI)<-"n0multiplierGrowthindividual"
	return(tI)
}

Migrationindividual<-function(collapseMatrix,complete=FALSE,n0multiplierMap,growthMap,migrationArray) {
	tI<-list(collapseMatrix=collapseMatrix,complete=complete,n0multiplierMap=n0multiplierMap,
		growthMap=growthMap,migrationArray=migrationArray)
	class(tI)<-"migrationindividual"
	return(tI)
}

CreateMSstringSpecific<-function(popVector,migrationIndividual,parameterVector,addedEventTime=NULL,
	addedEventTimeAsScalar=TRUE,nTrees=1,createSeed=TRUE){
	collapseMatrix<-migrationIndividual$collapseMatrix
	complete<-migrationIndividual$complete
	n0multiplierMap<-migrationIndividual$n0multiplierMap
	growthMap<-migrationIndividual$growthMap
	migrationArray<-migrationIndividual$migrationArray
	nPop<-length(popVector)
	numSteps<-dim(collapseMatrix)[2]
	n0multiplierParameters<-parameterVector[grep("n0multiplier",names(parameterVector))]
	growthParameters<-parameterVector[grep("growth",names(parameterVector))]
	migrationParameters<-parameterVector[grep("migration",names(parameterVector))]
	collapseParameters<-parameterVector[grep("collapse",names(parameterVector))]
	collapseCount<-1

	if(1>2 ){
	# stop(paste("Incorrect parameterVector: you passed ",length(parameterVector),"entries but it needs",
	#	length(parameterList),"entries:",paste(parameterList,sep=" ",collapse="")))
	}else{
		msString<-paste("-T -I",nPop,paste(popVector,sep=" ", collapse=" "),sep=" ")

		#do values at present
		initialN0multipler<-""
		initialGrowth<-""
		for (i in 1:nPop) {
			initialN0multipler<-paste(initialN0multipler,"-n",i,sprintf("%f",n0multiplierParameters[n0multiplierMap[i,1] ]),sep=" ")
			if(growthMap[i,1] > 0){
				initialGrowth<-paste(initialGrowth,"-g",i,sprintf("%f",growthParameters[growthMap[i,1] ]),sep=" ")
			}else{
				initialGrowth<-paste(initialGrowth,"-g",i,sprintf("%f",0),sep=" ")
			}
		}

		initialMigration<-" -ma"
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
		msString<-paste(msString,initialN0multipler,initialGrowth,initialMigration,sep="")

		#is there a collapse in the first gen?
		if(sum(!is.na(collapseMatrix[,1])) > 0){ #If this is not an added non-coalescence event
			if(max(collapseMatrix[,1],na.rm=TRUE) > 0) { #If there is a collapse
				initialCollapse<-""
				popsToCollapse<-which(collapseMatrix[,1]==1)
				for (i in 2:length(popsToCollapse)) {
					initialCollapse<-paste(initialCollapse, "-ej", sprintf("%f",collapseParameters[1]), popsToCollapse[i], popsToCollapse[1]  ,sep=" ")
				}
				msString<-paste(msString,initialCollapse,sep="")
			}
		}

		#now go back in time
		if (numSteps>1){
			eventCount<-1
			for (generation in 2:numSteps) {

				#If the previous generation was an added non-coalescence event with a hard-coded time period
				if(sum(!is.na(collapseMatrix[,generation - 1])) == 0){
					if(is.null(addedEventTime)){
						return(warning("Error: Need to specify an addedEventTime"))
					}
					if(addedEventTimeAsScalar == TRUE){
						if(sum(!is.na(collapseMatrix[,1:(generation - 1)])) > 0){ #If all the previous generations were NOT made up of added non-coalescence events
							collapseTime<-collapseParameters[collapseCount + 1] * addedEventTime[eventCount]
						}else{
							collapseTime<-collapseParameters[collapseCount] * addedEventTime[eventCount]
						}
					}else{
						collapseTime<-addedEventTime[eventCount]
					}
					eventCount<-eventCount + 1

				#Else use the values within the collapseParameter vector
				}else{
					collapseTime<-collapseParameters[collapseCount]
				}
				n0multiplierPositions<-which(!is.na(n0multiplierMap[,generation]))
				growthPositions<-which(!is.na(growthMap[,generation]))
				for (n0multiplierPosIndex in 1:length(n0multiplierPositions)) {
					n0multiplierPos<-n0multiplierPositions[n0multiplierPosIndex]
					msString<-paste(msString,"-en",sprintf("%f",collapseTime),sprintf("%f",n0multiplierPos),sprintf("%f",n0multiplierParameters[n0multiplierMap[n0multiplierPos,generation] ]),sep=" ")
				}
				for (growthPosIndex in 1:length(growthPositions)) {
					growthPos<-growthPositions[growthPosIndex]
					if(growthMap[growthPos,generation] > 0){
						msString<-paste(msString,"-eg",sprintf("%f",collapseTime),sprintf("%f",growthPos),sprintf("%f",growthParameters[growthMap[growthPos,generation] ]),sep=" ")
					}else{
						msString<-paste(msString,"-eg",sprintf("%f",collapseTime),sprintf("%f",growthPos),sprintf("%f",0),sep=" ")
					}
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
				if(sum(!is.na(collapseMatrix[,generation])) > 0){ #If the current generation is NOT an added non-coalescence event
					if(max(collapseMatrix[,generation],na.rm=TRUE) > 0){ #If the current generation contains a collapse event
						if(sum(!is.na(collapseMatrix[,1:(generation - 1)])) > 0){ #If all the previous generations were NOT made up of added non-coalescence events
				 			collapseCount<-collapseCount + 1 #then go to the next collapse parameter
						}
						initialCollapse<-""
						popsToCollapse<-which(collapseMatrix[,generation]>=1)
						for (i in 2:length(popsToCollapse)) {
							initialCollapse<-paste(initialCollapse, "-ej", sprintf("%f",collapseParameters[collapseCount]), popsToCollapse[i], popsToCollapse[1]  ,sep=" ")
						}
						msString<-paste(msString,initialCollapse,sep="")
					}
				}

			}
		}
		if (createSeed==TRUE){
			seed<-round(runif(3,1,10000000000))
			msString<-paste(msString,"-seed",seed[1],seed[2],seed[3],sep=" ")
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
		else if (ColCountLastZero(popIntervalsList[[i]]$collapseMatrix,0)==0) {
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
			if (length(which(abs(lastGen)>=1))>0) {
				survPop=survPop+1
			}else{
				stop("this population was complete")
			}
			rawIntervals<-c()
			if (survPop>1) {
				rawIntervals<-blockparts(c(1:survPop),survPop,include.fewer=TRUE)
				rawIntervals<-rawIntervals[,ColMax(rawIntervals)<=1] #we are okay with having all zeros: no population collapse
				rawIntervals<-rawIntervals[,ColCountIf(rawIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
			}else{
				rawIntervals=matrix(c(0,1),nrow=1)
			}
			rowMapping<-c(min(which(abs(lastGen)>=1)),which(lastGen==0))
			rawIntervalsRescaled<-matrix(NA,ncol=dim(rawIntervals)[2],nrow=length(lastGen))
			for (j in 1:length(rowMapping)) {
				rawIntervalsRescaled[rowMapping[j], ]<-rawIntervals[j,]
			}
			for (k in 1:dim(rawIntervalsRescaled)[2]){
				collapseValue<-max(currentParent$collapseMatrix[,dim(currentParent$collapseMatrix)[2]],na.rm=TRUE)
				rawIntervalsRescaled[,k][which(rawIntervalsRescaled[,k] == 1)]<-collapseValue + 1
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


#The basic idea here is that at each population in each time interval there is a n0multiplier. These can all be set to the same value,
#allowed to vary, or assigned in clumps (i.e., pops 1, 3, and 6 have the same n0multiplier value). This generates all such mappings,
#subject to staying within the maximum number of free parameters.
#Note that unlike for migration or growth, the minimum maximum number of n0multiplier parameters is one (instead of zero).
#This function now also adds in all possible growth models, which are generated using GenerateGrowthIndividuals
#The option has also been added to fix the n0muliplierMap or growthMap (or both) to be a particular matrix. These can be
#specified in the form of lists as n0multiplierList and growthList. Note that this must only be done when a single popInterval
#is specified (i.e., the collapse topology must also be fixed).
GenerateN0multiplierIndividuals<-function(popVector,popIntervalsList=GenerateIntervals(popVector),maxK=SetMaxK(popVector),maxN0K=Inf,
	maxGrowthK=Inf,n0multiplierList=NULL,growthList=NULL){
	n0multiplierIndividualsList<-list()
	for (i in sequence(length(popIntervalsList))){

		#If a particular n0multiplierMap is specified (in addition to a specified collapseMatrix), then just append this to each popInterval
		if(!is.null(n0multiplierList)){
			n0multiplierVec<-c()
			for(j in 1:length(n0multiplierList)){
				n0multiplierVec<-append(n0multiplierVec,n0multiplierList[[j]])
			}
			n0multiplierMap<-array(n0multiplierVec,dim=c(length(n0multiplierList[[1]]),length(n0multiplierList)))

			#Add in all possible growth models to each n0multiplierIndividual
			#If a particular growthMap is specified (in addition to a specified collapseMatrix), append it to each popInterval
			if(!is.null(growthList)){
				growthVec<-c()
				for(k in 1:length(growthList)){
					growthVec<-append(growthVec,growthList[[k]])
				}
				growthMap<-array(growthVec,dim=c(length(growthList[[1]]),length(growthList)))
				n0multiplierIndividualsList[[length(n0multiplierIndividualsList)+1]]<-N0multiplierGrowthindividual(popIntervalsList[[i]]$collapseMatrix,
					popIntervalsList[[i]]$complete,n0multiplierMap,growthMap)

			}else{

				#Else, do all possible growth maps
				growthMapList<-GenerateGrowthIndividuals(popVector=popVector,popInterval=popIntervalsList[[i]],maxK=maxK,maxGrowthK=maxGrowthK)
				for(k in sequence(length(growthMapList))){
					growthMap<-growthMapList[[k]]
					n0multiplierIndividualsList[[length(n0multiplierIndividualsList)+1]]<-N0multiplierGrowthindividual(popIntervalsList[[i]]$collapseMatrix,
						popIntervalsList[[i]]$complete,n0multiplierMap,growthMap)
				}
			}

		#Else, do all possible n0multiplier maps
		}else{
			n0multiplierMapTemplate<-1+0*popIntervalsList[[i]]$collapseMatrix  #will have all the populations, all with either NA or 1
			numLineages=sum(n0multiplierMapTemplate,na.rm=TRUE)
			possibleMappings<-AllParamCombinations(numCells=numLineages,minVal=1,maxVal=max(1,min(maxN0K,maxK-KPopInterval(popIntervalsList[[i]]))),
				allowableMaxInitial=1)
			for (mappingIndex in sequence(dim(possibleMappings)[1])) {
				thisMapping<-possibleMappings[mappingIndex,]
				n0multiplierMap<-n0multiplierMapTemplate
				whichPositions <- which(n0multiplierMap==1)
				n0multiplierMap[whichPositions]<-thisMapping

				#Add in all possible growth models to each n0multiplierIndividual
				#If a particular growthMap is specified (in addition to a specified collapseMatrix), append it to each popInterval
				if(!is.null(growthList)){
					growthVec<-c()
					for(k in 1:length(growthList)){
						growthVec<-append(growthVec,growthList[[k]])
					}
					growthMap<-array(growthVec,dim=c(length(growthList[[1]]),length(growthList)))
					n0multiplierIndividualsList[[length(n0multiplierIndividualsList)+1]]<-N0multiplierGrowthindividual(popIntervalsList[[i]]$collapseMatrix,
						popIntervalsList[[i]]$complete,n0multiplierMap,growthMap)

				}else{

					#Else, do all possible growth maps
					growthMapList<-GenerateGrowthIndividuals(popVector=popVector,popInterval=popIntervalsList[[i]],maxK=maxK,maxGrowthK=maxGrowthK)
					for(k in sequence(length(growthMapList))){
						growthMap<-growthMapList[[k]]
						n0multiplierIndividualsList[[length(n0multiplierIndividualsList)+1]]<-N0multiplierGrowthindividual(popIntervalsList[[i]]$collapseMatrix,
							popIntervalsList[[i]]$complete,n0multiplierMap,growthMap)
					}
				}
			}
		}
	}
	return(n0multiplierIndividualsList)
}

#For a given model including collapse and n0multiplier, this produces a list of all possible growth models
#accounting for the specified number of possible and maximum parameters
GenerateGrowthIndividuals<-function(popVector,popInterval,maxK=SetMaxK(popVector),maxGrowthK=Inf){
	growthIndividualsList<-list()
	growthMapTemplate<-1+0*popInterval$collapseMatrix  #will have all the populations, all with either NA or 1
	if(maxGrowthK == 0){
		growthMap<-growthMapTemplate * 0
		return(list(growthMap))
	}else{
		numLineages=sum(growthMapTemplate,na.rm=TRUE)
		possibleMappings<-AllParamCombinations(numCells=numLineages,minVal=0,maxVal=max(1,min(maxGrowthK,maxK-KPopInterval(popInterval))),
			allowableMaxInitial=1)
		for (mappingIndex in sequence(dim(possibleMappings)[1])) {
			thisMapping<-possibleMappings[mappingIndex,]
			growthMap<-growthMapTemplate
			whichPositions <- which(growthMap==1)
			growthMap[whichPositions]<-thisMapping
			growthIndividualsList[[length(growthIndividualsList)+1]]<-growthMap
		}
		return(growthIndividualsList)
	}
}

GetNumTreeHistories<-function(popVector,maxK=SetMaxK(popVector),maxN0K=1){
	maxK<-maxK + 1 #to account for fact that first n0 multiplier is not free
	popIntervalsList<-GenerateIntervals(popVector)
	n0multiplierIndividualsList<-GenerateN0multiplierIndividuals(popVector,popIntervalsList,maxK=maxK, maxN0K=maxN0K)
	return(length(n0multiplierIndividualsList))
}

#Now we will generate all possible assignments of pairwise migration. Again, we want to keep the total number of free parameters
#(times, n0multipliers, growth parameters, and migration rates) under our chosen max.
#Do we allow a model where migrations change anywhere along branch, or only at coalescent nodes?
#The problem with the latter is that you cannot fit some reasonable models: i.e., two populations persisting through time.
#Problem with the former is parameter space.

#In general, this function produces a list of all possible demographic models given a number of populations (specified by popVector),
#specified numbers of free parameters, and other optional filtering criteria. The list of models is call a migrationArray (a name which,
#confusingly, is also given to migration matrices within a model) and each model is referred to as a migrationIndividual.

#Each migrationIndividual is composed of four major demographic components:
#1. The coalescence history (the collapseMatrix, in which rows correspond to populations and columns correspond to temporal events)
#2. Population size scalar parameters (the n0multiplierMap, with the same dimensions as CollapseMatrix)
#3. Growth (alpha) parameters (the growthMap, also with the same dimensions as CollapseMatrix)
#4. Migration rate matrices (the migrationArray, where the number of matrices is equal to the number of temporal events specified by
#	the collapseMatrix.

#One can specify the total maximum number of free parameters (maxK), and also specify the maximum number of free parameters for a given
#parameter type (e.g., maxN0K, maxGrowthK, and maxMigrationK). Note that maxGrowthK and maxMigrationK can be set to zero; however the
#lowest possible value for N0K is one (in which case all population sizes are the same).

#Specific demographic components can also be fixed, such that the migrationArray generated only varies subsets of parameters. Fixed
#components are specified by collapseList, n0multiplierList, growthList, and migrationList. Note that a collapseList must always be
#specified in order to specify any other demographic component. The format for specifying collapseList, n0multiplierList, and
#growthList are the same: a list of vectors. A vector contains parameter values for each population (e.g., becoming the rows in
#collapseMatrix), and each vector in the list represents a different temporal event (e.g., becoming the columns in collapseMatrix).
#For example, collapseList = list(c(1,1,0),c(2,NA,2)) means that there are two coalescent events: in the first event, population 1 and 2
#coalesce while population 3 does not; in the second event, ancestral population 1-2 coalesces with population 3.

#MigrationLists differ in that a list of matrices, rather than vectors, must be specified. There will be one migration matrix
#for the present generation, plus any historical matrices that apply (there will be as many matrices as there are collapse events).
#So, in the three population scenario, if there is symmetrical migration between populations 1 and 2 and no historical migration,
#migrationList will be:

# migrationList<-list(
# t(array(c(
# NA, 1, 0,
# 1, NA, 0,
# 0, 0, NA),
# dim=c(3,3))),
#
# t(array(c(
# NA, NA, 0,
# NA, NA, NA,
# 0,  NA, NA),
# dim=c(3,3))))

#Note that in R, arrays are constructed by reading in values from column 1 first then from column 2, etc. However, it is more intuitive
#to construct migration matrices by rows (i.e., first listing values for row 1, then row 2, etc). Thus, I think it is easier to type in the
#arrays by rows (as done above), and then transpose them (using "t"). Also, spacing and hard returns can be used to visualize these values
#in the form of matrices.

#Other methods of model filtering (using forceSymmetricalMigration or forceTree) can also be implemented.
GenerateMigrationIndividuals<-function(popVector,maxK=SetMaxK(popVector),maxN0K=1,maxGrowthK=0,maxMigrationK=1,
	collapseList=NULL,n0multiplierList=NULL,growthList=NULL,migrationList=NULL,forceSymmetricalMigration=TRUE,
	forceTree=FALSE,verbose=FALSE,parallelRep=NULL){

	#Check for compatibilities among inputs
	if((!is.null(n0multiplierList) || !is.null(growthList) || !is.null(migrationList)) && is.null(collapseList)){
		return(message("Error: when specifying n0multiplierList, growthList, or migrationLists, a collapseList must also be specified."))
	}

	#Some preliminaries
	maxK<-maxK + 1 #to account for fact that first n0 multiplier is not free
	migrationIndividualsList<-list()

	#If a fixed collapse scenario is specified, use that
	if(!is.null(collapseList)){
		collapseVec<-c()
		for(j in 1:length(collapseList)){
			collapseVec<-append(collapseVec,collapseList[[j]])
		}
		collapseMatrix<-array(collapseVec,dim=c(length(collapseList[[1]]),length(collapseList)))
		popIntervalsList<-list(Popinterval(collapseMatrix,complete=TRUE))
	}else{

	#Else, do all possible collapse scenarios
		popIntervalsList<-GenerateIntervals(popVector)
	}

	#Copy these models across various possible n0multiplier and growth parameter scenarios
	n0multiplierIndividualsList<-GenerateN0multiplierIndividuals(popVector=popVector,popIntervalsList=popIntervalsList,
		maxK=maxK,maxN0K=maxN0K,maxGrowthK=maxGrowthK,n0multiplierList=n0multiplierList,growthList=growthList)

	#Run for a specific tree history or for all at once?
	if(!is.null(parallelRep)){
		treeHistoryNum<-parallelRep
	}else{
		treeHistoryNum<-sequence(length(n0multiplierIndividualsList))
	}
	for (i in treeHistoryNum) {
		if (verbose==TRUE) {
			print(paste("doing ",i,"/",length(n0multiplierIndividualsList)))
		}
		collapseMatrix<-n0multiplierIndividualsList[[i]]$collapseMatrix
		n0multiplierMap<-n0multiplierIndividualsList[[i]]$n0multiplierMap
		growthMap<-n0multiplierIndividualsList[[i]]$growthMap
		numFinalPops<-dim(collapseMatrix)[1]
		numSteps<-dim(collapseMatrix)[2]

		#Check for input compatibilities
		if((paste(dim(collapseMatrix),collapse="") != paste(dim(n0multiplierMap),collapse="")) ||
		(paste(dim(collapseMatrix),collapse="") != paste(dim(growthMap),collapse="")) ||
		(paste(dim(growthMap),collapse="") != paste(dim(n0multiplierMap),collapse=""))){
			return(message("Error: specified collapseLists, n0multiplierLists, and growthLists must have the same dimensions"))
		}

		if((KCollapseMatrix(collapseMatrix) + KN0multiplierMap(n0multiplierMap) + KGrowthMap(growthMap) ) <= maxK){

			#If a specific migrationArray is specified, just use that one
			if(!is.null(migrationList)){
				migrationVec<-c()
				for(i in 1:length(migrationList)){
					migrationVec<-append(migrationVec,migrationList[[i]])
				}
				migrationArray<-array(migrationVec,dim=c(nrow(migrationList[[1]]),ncol(migrationList[[1]]),
					length(migrationList)))

				#Check for input compatibilities
				if(dim(collapseMatrix)[1] != dim(migrationArray)[1] || dim(collapseMatrix)[2] != dim(migrationArray)[3]){
					return(message("Error: specified collapseLists and migrationLists must have the same dimensions"))
				}

				newIndividual<-Migrationindividual(collapseMatrix,n0multiplierIndividualsList[[i]]$complete,
					n0multiplierMap,growthMap,migrationArray)
				#Toss models in which genes cannot coalesce or which exceed the max number of parameters
				if(CheckFiniteCoalescence(newIndividual) && KAll(newIndividual) < maxK){
					migrationIndividualsList[[length(migrationIndividualsList)+1]]<-newIndividual
				}

			#Otherwise...
			}else{

				#If no migration desired, just create empty matrices
				if(maxMigrationK == 0){
					migrationArray<-array(data=NA,dim=c(nrow(collapseMatrix),nrow(collapseMatrix),ncol(collapseMatrix))) #dimnames=c("from","to","generation")
					for(interval in 1:ncol(collapseMatrix)) {
						for(fromPop in 1:nrow(collapseMatrix)) {
							for(toPop in 1:nrow(collapseMatrix)) {
								if(fromPop!=toPop) {
									if(!is.na(collapseMatrix[fromPop,interval]) && !is.na(collapseMatrix[toPop,interval])) {
										migrationArray[fromPop,toPop,interval]<-0
									}
								}else{
									if(!is.na(collapseMatrix[fromPop,interval]) && !is.na(collapseMatrix[toPop,interval]) &&
										toPop>fromPop) { #so only fill top diagonal
										migrationArray[fromPop,toPop,interval]<-0
									}
								}
							}
						}
					}
					newIndividual<-Migrationindividual(collapseMatrix,n0multiplierIndividualsList[[i]]$complete,
						n0multiplierMap,growthMap,migrationArray)
					#Toss models in which genes cannot coalesce or which exceed the max number of parameters
					if(CheckFiniteCoalescence(newIndividual) && KAll(newIndividual) < maxK){
						migrationIndividualsList[[length(migrationIndividualsList)+1]]<-newIndividual
					}
				}else{

					#Else, produce all possible matrices
					migrationTemplate<-array(data=NA,dim=c(numFinalPops,numFinalPops,numSteps)) #dimnames=c("from","to","generation")
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
					mappingMaxKAllowed <- min(maxK - ( KCollapseMatrix(collapseMatrix) + KN0multiplierMap(n0multiplierMap) ),maxMigrationK)

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
					while(evaluateMapping){
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
						newIndividual<-Migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix,n0multiplierIndividualsList[[i]]$complete,
							n0multiplierMap,growthMap,migrationArray)
						if(CheckFiniteCoalescence(newIndividual) && KAll(newIndividual)<maxK){
							migrationIndividualsList[[length(migrationIndividualsList)+1]]<-newIndividual
						}
						#Code for storing models separately for each tree history
		#				if(CheckFiniteCoalescence(newIndividual) && KAll(newIndividual)<maxK  && parameterize==TRUE){
		#					if(!exists(paste("migrationIndividualsList_",i,sep=""))){
		#						assign(paste("migrationIndividualsList_",i,sep=""),list())
		#						#Method for adding elements to dynamically named list
		#						assign(paste("migrationIndividualsList_",i,sep=""),
		#						{x<-get(paste("migrationIndividualsList_",i,sep=""))
		#						x[[length(get(paste("migrationIndividualsList_",i,sep=""))) + 1]]<-4
		#						x})
		#					}
		#				}
						numevaluations<-numevaluations+1
						if(numevaluations<=dim(allMappings)[1]) {
							thisMapping<-allMappings[numevaluations,]
						} else {
							evaluateMapping<-FALSE
						}
					}#end while
				}#end else (zero versus non-zero parameters)
			}#end else (specified matrix versus non-specified matrix)
		}#end if
	}#end for
	if(forceTree) {
		migrationIndividualsList<-FilterForFullyResolved(migrationIndividualsList)
	}
	return(migrationIndividualsList)
}

# Previous versions of PHRAPL produced model lists (migrationArrays) that lack growth
# matrices (growthMaps). However, growthMaps must now be present to analyze models in PHRAPL,
# even if growth parameters are not incorporated in the model. Thus, to facilitate the analysis
# of old migrationArrays, this function takes as input a  migrationArray lacking growthMaps
# and adds an empty growthMap to each model (i.e., a growthMap filled with zeros) within the
# migrationArray.
AddGrowthToAMigrationArray<-function(migrationArray){
	if(class(migrationArray) == "migrationindividual"){
		migrationArray<-list(migrationArray)
	}
	for(i in 1:length(migrationArray)){
		growthMap<-migrationArray[[i]]$collapseMatrix
		growthMap[!is.na(growthMap)]<-0
		migrationArray[[i]]<-Migrationindividual(migrationArray[[i]]$collapseMatrix,
			migrationArray[[i]]$complete,migrationArray[[i]]$n0multiplierMap,
			growthMap,migrationArray[[i]]$migrationArray)
	}
	return(migrationArray)
}

# This function integrates a non-coalescence demographic event within a model
# or set of models. This can be useful if one wants to posit shifts in a parameter that
# do not correspond to a splitting event (e.g., one would like migration to only occur
# for a given length of time after a split, but then to cease).
#
# To use this function, one must specify a model (migrationIndividual) or set of models
# (migrationArray). If a set of models is specified, these models must all contain the
# same number of populations and the same number of collapse events (i.e., the collapseMatrix
# compenent within each migrationIndividual must have the same dimensions).
#
# The relative timing of the desired new event must be specified as a single number using
# eventTime. An eventTime of 1 will place a new event (i.e., column) in the first column
# position within the collapseMatrix (and the other columns will be shifted to the right).
# An eventTime of 2 will place the new column in the second position, etc. The eventTime
# cannot exceed the number of events (i.e., columns) within the original collapseMatrix.
# The added column will consist entierly of NAs, indicating that no population coalescence
# can occur at this time.
#
# Finally, one can specify a new set of n0multiplier, growth, and/or migration parameter
# indexes to be invoked at the new specified time period using n0muliplierVec, growthVec, or
# migrationMat, respectively. When the default value for these is used (which is NULL),
# n0multiplier, growth, and migration matrices within each model are automatically expanded
# by simply copying over parameter indexes from the adjacent time period to the new time
# period (i.e., no change is invoked). For n0multiplier and growth, a new column is added;
# for migration, a new matrix is added.
#
# However, one can also specify the parameter indexes to be used at this new time, allowing
# a shift in that parameter, without a correponding coalescence event. These can either be
# specified as a single value, which will be applied to all populations, or as vectors
# (for n0multiplier, growth, or migration) or matrices (for migration only) containing the
# relevant parameter indexes for each population (NAs can be included or excluded).
AddEventToMigrationArray<-function(migrationArray,eventTime,n0multiplierVec=NULL,growthVec=NULL,migrationMat=NULL){
	if(class(migrationArray) == "migrationindividual"){
		migrationArray<-list(migrationArray)
	}
	for(i in 1:length(migrationArray)){
		collapseMatrixOld<-migrationArray[[i]]$collapseMatrix
		n0multiplierMapOld<-migrationArray[[i]]$n0multiplierMap
		growthMapOld<-migrationArray[[i]]$growthMap
		migrationArrayOld<-migrationArray[[i]]$migrationArray

		if(eventTime > dim(collapseMatrixOld)[2]){
			return(warning("Error: The specified eventTime exceeds the current number of collapse events"))
		}

		#Insert new time event
		collapseMatrix<-matrix(NA,nrow=dim(collapseMatrixOld)[1],ncol=(dim(collapseMatrixOld)[2] + 1))
		adjust<-0
		for(colTime in 1:(dim(collapseMatrixOld)[2] + 1)){
			if(eventTime != colTime){
				collapseMatrix[,colTime]<-collapseMatrixOld[,(colTime - adjust)]
			}else{
				adjust<-1
			}
		}

		#Insert new n0multiplier
		n0multiplierMap<-matrix(NA,nrow=dim(n0multiplierMapOld)[1],ncol=(dim(n0multiplierMapOld)[2] + 1))
		n0multiplierMapOldReset<-n0multiplierMapOld
		adjust<-0
		for(n0Time in 1:(dim(n0multiplierMapOld)[2] + 1)){
			n0multiplierMapOld<-n0multiplierMapOldReset
			if(eventTime != n0Time){
				n0multiplierMap[,n0Time]<-n0multiplierMapOld[,(n0Time - adjust)]
			}else{
				if(!is.null(n0multiplierVec)){
					if(length(n0multiplierVec) == 1){
						n0multiplierVecNew<-n0multiplierMapOld[,n0Time]
						n0multiplierVecNew[!is.na(n0multiplierVecNew)]<-n0multiplierVec
						n0multiplierMap[,n0Time]<-n0multiplierVecNew
					}else{
						if(length(n0multiplierMapOld[,n0Time][!is.na(n0multiplierMapOld[,n0Time])]) !=
							length(n0multiplierVec[!is.na(n0multiplierVec)])){
							return(warning("Error: The number of n0multiplier values specified does not match the desired history"))
						}
						n0multiplierMapOld[,n0Time][!is.na(n0multiplierMapOld[,n0Time])]<-
							n0multiplierVec[!is.na(n0multiplierVec)]
						n0multiplierMap[,n0Time]<-n0multiplierMapOld[,n0Time]
					}
				}else{
					n0multiplierMap[,n0Time]<-n0multiplierMapOld[,n0Time]
				}
				adjust<-1
			}
		}

		#Insert new growth
		growthMap<-matrix(NA,nrow=dim(growthMapOld)[1],ncol=(dim(growthMapOld)[2] + 1))
		growthMapOldReset<-growthMapOld
		adjust<-0
		for(gTime in 1:(dim(growthMapOld)[2] + 1)){
			growthMapOld<-growthMapOldReset
			if(eventTime != gTime){
				growthMap[,gTime]<-growthMapOld[,(gTime - adjust)]
			}else{
				if(!is.null(growthVec)){
					if(length(growthVec) == 1){
						growthVecNew<-growthMapOld[,gTime]
						growthVecNew[!is.na(growthVecNew)]<-growthVec
						growthMap[,gTime]<-growthVecNew
					}else{
						if(length(growthMapOld[,gTime][!is.na(growthMapOld[,gTime])]) !=
							length(growthVec[!is.na(growthVec)])){
							return(warning("Error: The number of growth values specified does not match the desired history"))
						}
						growthMapOld[,gTime][!is.na(growthMapOld[,gTime])]<-
							growthVec[!is.na(growthVec)]
						growthMap[,gTime]<-growthMapOld[,gTime]
					}
				}else{
					growthMap[,gTime]<-growthMapOld[,gTime]
				}
				adjust<-1
			}
		}

		#Insert new migration
		numberOfValuesPerMatrix<-length(migrationArrayOld) / dim(migrationArrayOld)[3]
		migrationArrayNew<-array(rep(NA,(length(migrationArrayOld) + numberOfValuesPerMatrix)),
			dim=c(sqrt(numberOfValuesPerMatrix),sqrt(numberOfValuesPerMatrix),(dim(migrationArrayOld)[3] + 1)))
		migrationArrayOldReset<-migrationArrayOld
		adjust<-0

		for(migTime in 1:(dim(migrationArrayOld)[3] + 1)){
			migrationArrayOld<-migrationArrayOldReset
			if(eventTime != migTime){
				migrationArrayNew[,,migTime]<-migrationArrayOld[,,(migTime - adjust)]
			}else{
				if(!is.null(migrationMat)){
					if(length(migrationMat) == 1){
						migrationMatNew<-migrationArrayOld[,,migTime]
						migrationMatNew[!is.na(migrationMatNew)]<-migrationMat
						migrationArrayNew[,,migTime]<-migrationMatNew
					}else{
						if(length(migrationArrayOld[,,migTime][!is.na(migrationArrayOld[,,migTime])]) !=
							length(migrationMat[!is.na(migrationMat)])){
							return(warning("Error: The number of migration values specified does not match the desired history"))
						}
						migrationArrayOld[,,migTime][!is.na(migrationArrayOld[,,migTime])]<-migrationMat[!is.na(migrationMat)]
						migrationArrayNew[,,migTime]<-migrationArrayOld[,,migTime]
					}

				}else{
					migrationArrayNew[,,migTime]<-migrationArrayOld[,,migTime]
				}
				adjust<-1
			}
		}

		#Update migrationIndividual components
		migrationArray[[i]]$collapseMatrix<-collapseMatrix
		migrationArray[[i]]$n0multiplierMap<-n0multiplierMap
		migrationArray[[i]]$growthMap<-growthMap
		migrationArray[[i]]$migrationArray<-migrationArrayNew
	}
	return(migrationArray)
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
	if(!is.null(migrationIndividualsToKill)){
		migrationArray<-migrationArray[-1*migrationIndividualsToKill]
	}
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

GenerateMigrationIndividualsOneAtATime<-function(collapseList,n0multiplierList=NULL,growthList=NULL,migrationList=NULL){
	#The purpose of this function is create a single a priori migrationIndividual
	#for a given collapse, n0multiplier, growth, and migration history.
	#The four inputs for this function are "collapseList", "n0multiplierMap",
	#"growthMap", and "migrationList".

	#CollapseList is the only set of parameters that must be specified. This gives a
	#list of collapse history vectors, one vector for each
	#coalescent event in the tree. So, collapseList = list(c(1,1,0),c(2,NA,2)) means
	#that there are two coalescent events: in the first event, population 1 and 2
	#coalesce while population 3 does not; in the second event, ancestral population
	#1-2 coalesces with population 3.

	#The remaining three parameters may be specified or not specified. If not specified,
	#(i.e., set to NULL),then null matrices will be automatically constructed in which
	#all n0multipliers are set to one, and all growth and migration parameters are set to zero.

	#If specifying n0multiplierList and/or growthList, the format for these is the same as for
	#collapseList, and the available parameters for these must match the splitting history
	#depicted in the collapseList.

	#MigrationList is a list of migration matrices. There will be one migration matrix
	#for the present generation, plus any historical matrices that apply
	#(there will be as many matrices as there are collapse events).
	#So, in the three population scenario, if there is symmetrical migration between
	#populations 1 and 2 and no historical migration, migrationList will be:

	# migrationList<-list(
	# t(array(c(
	# NA, 1, 0,
	# 1, NA, 0,
	# 0, 0, NA),
	# dim=c(3,3))),
	#
	# t(array(c(
	# NA, NA, 0,
	# NA, NA, NA,
	# 0,  NA, NA),
	# dim=c(3,3))))

	#Note that in R, arrays are constructed by reading in values from column 1 first then from
	#column 2, etc. However, it is more intuitive to construct migration matrices by rows (i.e.,
	#first listing values for row 1, then row 2, etc). Thus, I think it is easier to type in the
	#arrays by rows (as done above), and then transpose them (using "t"). Also, spacing and hard
	#returns can be used to visualize these values in the form of matrices.

	#Create collapase matrix from list
	collapseVec<-c()
	for(i in 1:length(collapseList)){
		collapseVec<-append(collapseVec,collapseList[[i]])
	}
	collapseMatrix<-array(collapseVec,dim=c(length(collapseList[[1]]),length(collapseList)))

	#Make n0multiplierMap
	if(is.null(n0multiplierList)){
		collapseVec[!is.na(collapseVec)]<-1
		n0multiplierMap<-array(collapseVec,dim=c(length(collapseList[[1]]),length(collapseList)))
	}else{
		n0multiplierVec<-c()
		for(i in 1:length(n0multiplierList)){
			n0multiplierVec<-append(n0multiplierVec,n0multiplierList[[i]])
		}
		n0multiplierMap<-array(n0multiplierVec,dim=c(length(n0multiplierList[[1]]),length(n0multiplierList)))
	}

	#Make growthMap
	if(is.null(growthList)){
		collapseVec[!is.na(collapseVec)]<-0
		growthMap<-array(collapseVec,dim=c(length(collapseList[[1]]),length(collapseList)))
	}else{
		growthVec<-c()
		for(i in 1:length(growthList)){
			growthVec<-append(growthVec,growthList[[i]])
		}
		growthMap<-array(growthVec,dim=c(length(growthList[[1]]),length(growthList)))
	}

	#Make migrationArray
	if(is.null(migrationList)){
		migrationArray<-array(data=NA,dim=c(nrow(collapseMatrix),nrow(collapseMatrix),ncol(collapseMatrix)))
		for(interval in 1:ncol(collapseMatrix)) {
			for(fromPop in 1:nrow(collapseMatrix)) {
				for(toPop in 1:nrow(collapseMatrix)) {
					if(fromPop!=toPop) {
						if(!is.na(collapseMatrix[fromPop,interval]) && !is.na(collapseMatrix[toPop,interval])) {
							migrationArray[fromPop,toPop,interval]<-0
						}
					}else{
						if(!is.na(collapseMatrix[fromPop,interval]) && !is.na(collapseMatrix[toPop,interval]) &&
							toPop>fromPop) { #so only fill top diagonal
							migrationArray[fromPop,toPop,interval]<-0
						}
					}
				}
			}
		}
	}else{
		migrationVec<-c()
		for(i in 1:length(migrationList)){
			migrationVec<-append(migrationVec,migrationList[[i]])
		}
		migrationArray<-array(migrationVec,dim=c(nrow(migrationList[[1]]),ncol(migrationList[[1]]),
			length(migrationList)))
	}

	migrationIndividual<-Migrationindividual(collapseMatrix=collapseMatrix,complete=TRUE,
		n0multiplierMap=n0multiplierMap,growthMap=growthMap,migrationArray=migrationArray)
	return(migrationIndividual)
}

LoadMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation=system.file("msdir","ms",package="P2C2M")) {
	msCallInfo<-CreateMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
	geneTrees<-system(paste(msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';'"),intern=TRUE)
	return(geneTrees)
}

SaveMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation=system.file("msdir","ms",package="P2C2M"),file="sim.tre") {
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

PrepSubsampling<-function(assignmentsGlobal,observedTrees,popAssignments,subsamplesPerGene,nIndividualsDesired=NULL,minPerPop=1,
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
				warning(paste("Warning: Tree number ",tree,"contains tip names not included in the inputted assignment file.",
					"These tips will not be subsampled.", sep=" "))
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
	subsampleWeightsList<-list()
	numTrees<-0
	treeCounter<-0
	#Get number of observed trees
	for(t in sequence(length(observedTrees))){
		if (class(observedTrees[[t]]) != "multiPhylo"){
        	observedTrees[[t]] <- c(observedTrees[[t]])
    	}
		numTrees<-numTrees + length(observedTrees[[t]])
	}
	for(f in sequence(length(observedTrees))){ #for each popAssignments
	subsampleWeights<-data.frame(weight=matrix(NA,ncol=1,nrow=length(observedTrees[[f]])))
		for(y in sequence(length(observedTrees[[f]]))){ #for each tree
			subsampleWeights[y,1]<-GetPermutationWeights(phy=observedTrees[[f]][[y]],
				popVector=popAssignments[[f]])
			treeCounter<-treeCounter + 1
			cat(paste("Finished with ",treeCounter," out of ",numTrees," trees\n",sep=""))
		}
	subsampleWeightsList[[f]]<-subsampleWeights
	}
	return(subsampleWeightsList)
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
	migrationArrayMap<-data.frame(matrix(c(1,1,1,1,1,1,1,1,1),nrow=1)) #uses the first entry in migrationArray
	for (model in 2:length(migrationArray)) {
		newRow<-c(model,rep(NA,8))
    	for (comparison in sequence(dim(migrationArrayMap)[1])) {
      		comparisonRow<-as.integer(migrationArrayMap[comparison,])
     		if(identical(migrationArray[[comparisonRow[3] ]]$collapseMatrix, migrationArray[[model]]$collapseMatrix)) {
       			newRow[2:3]<-comparisonRow[2:3]
     		}
    		if(identical(migrationArray[[comparisonRow[5] ]]$n0multiplierMap, migrationArray[[model]]$n0multiplierMap)) {
       			newRow[4:5]<-comparisonRow[4:5]
    		}
     		if(identical(migrationArray[[comparisonRow[7] ]]$growthMap, migrationArray[[model]]$growthMap)) {
       			newRow[4:5]<-comparisonRow[6:7]
     		}
     		if(identical(migrationArray[[comparisonRow[9] ]]$migrationArray, migrationArray[[model]]$migrationArray)) {
       			newRow[6:7]<-comparisonRow[8:9]
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
		if(is.na(newRow[8])) {
      		newRow[8:9]<-c(1+max(migrationArrayMap[,8]),model)
    	}
    	migrationArrayMap<-rbind(migrationArrayMap,newRow)
  	}
  names(migrationArrayMap)<-c("model","collapseMatrix.number","collapseMatrix.parent","n0multiplierMap.number","n0multiplierMap.parent",
  	"growthMap.number","growthMap.parent","migrationArray.number","migrationArray.parent")
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
#of models into a table (nloptr or genoud)
ExtractOptimizationAICs<-function(result=result,migrationArray=migrationArray,modelRange=c(1:length(migrationArray)),
	setCollapseZero=NULL,optimization=NULL){

	#Pull out the overall AICs
	if(optimization == "nloptr"){
		AIC<-grep("objective",result,value=TRUE)
		AIC<-gsub(".+objective = ","",AIC)
		AIC<-gsub(", solution.+","",AIC)
	}
	if(optimization == "genoud"){
		AIC<-grep("value",result,value=TRUE)
		AIC<-gsub(".+value = ","",AIC)
		AIC<-gsub(", par.+","",AIC)
	}

	#Construct dataframe consisting of AICs and model descriptions for each model
  	overall.results<-data.frame
	params.K<-rep(NA, length(AIC))
	params.vector<-rep(NA, length(AIC))
	for (i in sequence(length(AIC))) {
		#for K, subtract off fixed values
		params.K[i]<-(length(MsIndividualParameters(migrationArray[[i]])) - 1) - length(setCollapseZero)
		params.vector[i]<-paste(MsIndividualParameters(migrationArray[[i]]), sep=" ", collapse=" ")
		if(!is.null(setCollapseZero)){ #toss collapse parameter names when they are set to zero
			params.vector[i]<-paste(strsplit(params.vector[i]," ")[[1]][-setCollapseZero],sep=" ",collapse=" ")
		}
		#Remove first n0multiplier
        if(length(grep("n0multiplier_1",params.vector[i])) > 0){
        	params.vector[i] <- paste(strsplit(params.vector[i]," ")[[1]][-(grep("n0multiplier",
        		strsplit(params.vector[i]," ")[[1]])[1])], seo = " ", collapse = "")
        }
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
#of models into a table (grid)
ExtractGridAICs<-function(result=result,migrationArray=migrationArray,modelRange=c(1:length(migrationArray)),
	setCollapseZero=NULL){

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
		params.vector[i]<-paste(colnames(result[[i]])[-1], sep=" ", collapse=" ")
		#for K, subtract off fixed values
		if(length(grep("n0multiplier_1",params.vector[i])) > 0){
			params.K[i]<-(ncol(result[[i]]) - 1) - length(setCollapseZero) - 1
		}else{
			params.K[i]<-(ncol(result[[i]]) - 1) - length(setCollapseZero)
		}
		if(!is.null(setCollapseZero)){ #toss collapse parameter names when they are set to zero
			params.vector[i]<-paste(strsplit(params.vector[i]," ")[[1]][-setCollapseZero],sep=" ",collapse=" ")
		}
		#Remove first n0multiplier
        if(length(grep("n0multiplier_1",params.vector[i])) > 0){
        	params.vector[i] <- paste(strsplit(params.vector[i]," ")[[1]][-(grep("n0multiplier",
        		strsplit(params.vector[i]," ")[[1]])[1])], seo = " ", collapse = "")
        }
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







###############################Post-analysis Functions####################


#This function takes output from an exhaustive search (using nLoptr optimization), and assembles parameter
#indexes and estimates based on a set of models. A list containing two tables is outputted: one containing
#parameter indexesand one containing parameter estimates
#Note that this function produces some parameters that entail ambiguous names (for example, for ancestral
#migration rates, it can't distinguish different topologies). Thus, this method of calculating parameters
#has been abandoned (December 2015) in favor of an unambiguous naming method, implemented using the
#function ExtractUnambiguousNLoptrParameters. However, this old method can be implemented in a GridSearch
#by setting parameter_ambiguous=FALSE.
ExtractParameters<-function(migrationArray=migrationArray,result=result,popVector,rm.n0=TRUE){

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
			#Now take away from the parameters vector the used up parameter values
			currentSol<-tail(currentSol,length(currentSol) - uniqueN0multi)
		}




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
	if(rm.n0==TRUE){
		parameterIndexes<-cbind(collapseParms,migParms)
	}else{
		parameterIndexes<-cbind(collapseParms,n0multiParms,migParms)
	}
	#Add an "I" to colnames for the indexes (to distinguish them from the parameter values colnames)
	for(rep in 1:ncol(parameterIndexes)){
		colnames(parameterIndexes)[rep]<-paste(colnames(parameterIndexes)[rep],"_I",sep="")
	}

	if(rm.n0==TRUE){
		parameterValues<-cbind(collapseParmsValues,migParmsValues)
	}else{
		parameterValues<-cbind(collapseParmsValues,n0multiParmsValues,migParmsValues)
	}
	parameterAll<-cbind(parameterValues,parameterIndexes)

	return(list(parameterValues,parameterIndexes))
}

#This function takes output from a grid search and assembles parameter
#indexes and estimates based on a set of models. A list containing two tables is outputted: one containing
#parameter indexesand one containing parameter estimates
#Note that this function produces some parameters that entail ambiguous names (for example, for ancestral
#migration rates, it can't distinguish different topologies). Thus, this method of calculating parameters
#has been abandoned (December 2015) in favor of an unambiguous naming method, implemented using the
#function ExtractUnambiguousGridParameters. However, this old method can be implemented in a GridSearch
#by setting parameter_ambiguous=FALSE.
ExtractGridParameters<-function(migrationArray=migrationArray,result=result,popVector,rm.n0=TRUE,dAIC.cutoff=2){

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


	#############FILL IN A MATRIX WITH THE COLLAPSE  AND N0MULTI PARAMETER INDEXES

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
		#Make vector of top parameter estimates from the grid (based on the mean of the top five estimates for each parm)
		currentSol<-c()
		positionOfNonParameters <- grep("lnL", colnames(result[[eachCol]]))[1] #find slopes
		if(!is.na(positionOfNonParameters)){
			currentResult<-result[[eachCol]][2:(positionOfNonParameters - 1)] #toss slopes
		}else{
			currentResult<-result[[eachCol]][2:length(result[[eachCol]])] #if only 1 popAssignment used
		}

		#For each column of parameters, get optimal parameter values according the define dAIC cutoff
		for(rep in 1:ncol(currentResult)){
			dAIC.temp<-c()
			for(a in 1:nrow(result[[eachCol]][1])){
				dAIC.temp<-append(dAIC.temp,result[[eachCol]][a,1] - result[[eachCol]][1,1])
			}
			currentSol<-append(currentSol,mean(currentResult[which(dAIC.temp < dAIC.cutoff),rep])) #Note that grid prams are already back-transformed (i.e., logged)
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
			#Now take away from the parameters vector the used up parameter values
			currentSol<-tail(currentSol,length(currentSol) - uniqueN0multi)
		}


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
	if(rm.n0==TRUE){
		parameterIndexes<-cbind(collapseParms,migParms)
	}else{
		parameterIndexes<-cbind(collapseParms,n0multiParms,migParms)
	}
	#Add an "I" to colnames for the indexes (to distinguish them from the parameter values colnames)
	for(rep in 1:ncol(parameterIndexes)){
		colnames(parameterIndexes)[rep]<-paste(colnames(parameterIndexes)[rep],"_I",sep="")
	}

	if(rm.n0==TRUE){
		parameterValues<-cbind(collapseParmsValues,migParmsValues)
	}else{
		parameterValues<-cbind(collapseParmsValues,n0multiParmsValues,migParmsValues)
	}
	parameterAll<-cbind(parameterValues,parameterIndexes)

	return(list(parameterValues,parameterIndexes))
}

#For a given phrapl result object and migrationArray, this function makes a table of
#model averaged parameter values for each parameter and model. The parameters are
#unambiguous in that tip population information is included in ancestral population names.
#This function is implemented at the end of a GridSearch and is called when parameters
#are estimated using a parameter grid.
ExtractUnambiguousGridParameters<-function(overall.results=NULL,gridList=NULL,migrationArray,
	rm.n0=TRUE,longNames=FALSE,sortParameters=TRUE,sortModelsAIC=TRUE,nonparmCols=2){

	#Add parameters
	all.parameters <- unique(sort(sapply(migrationArray,
		MsIndividualParametersConversionToUnambiguous,unambiguous.only = TRUE)))
	results.temp <- c()
	row.count = 0
	#Get models
	whichGrid<-1
	for (models in overall.results$models){ #Get current models
		param.mappings <- MsIndividualParametersConversionToUnambiguous(migrationArray[[whichGrid]],
			longNames=longNames)
		results.temp <- rbind(results.temp, t(apply(gridList[[whichGrid]],
			MARGIN = 1, DoSingleRowExtraction, param.mappings = param.mappings,
			all.parameters = all.parameters, models = models)))
	   whichGrid<-whichGrid + 1
	}
	#Convert all zero values to NA
	#Note that if tau = NA, this means tau was set to zero in the analysis and thus should stay
	#that way for model averaging
	results.temp<-data.frame(results.temp, check.names=FALSE, row.names=c(1:nrow(results.temp)))
	if(longNames==TRUE){
		results.temp[, !grepl("collapse",colnames(results.temp))][,-(0:nonparmCols)][results.temp[,
			!grepl("collapse",colnames(results.temp))][,-(0:nonparmCols)] == 0] <- NA
	}else{
		results.temp[, !grepl("t",colnames(results.temp))][,-(0:nonparmCols)][results.temp[,
			!grepl("t",colnames(results.temp))][,-(0:nonparmCols)] == 0] <- NA
	}
	#Remove undesired parameters
	if(rm.n0==TRUE){
		if(longNames==TRUE){
			if(length(grep("0multiplier", colnames(results.temp))) <= 1){
    			results.temp<-results.temp[,!grepl("n0multiplier", colnames(results.temp))]
    		}
    	}else{
    		if(length(grep("n", colnames(results.temp))) <= 1){
    			results.temp<-results.temp[,!grepl("n", colnames(results.temp))]
    		}
    	}
	}

	#Model average parameters
	modelAveragedResults<-ModelAverageEachModel(totalData=results.temp,parmStartCol=(nonparmCols + 1))

	#Sort parameter columns into desired order
	if(sortParameters == TRUE){
		modelAveragedResults<-sortParameterTable(modelAveragedResults,(nonparmCols + 1),longNames=longNames)
	}
    #Sort models
    if(sortModelsAIC == TRUE){
    	modelAveragedResults<-modelAveragedResults[order(modelAveragedResults$AIC),]
    }
    return(modelAveragedResults)
}

#For a given phrapl result object and migrationArray, this function makes a table of
#model averaged parameter values for each parameter and model. The parameters are
#unambiguous in that tip population information is included in ancestral population names.
#This function is implemented at the end of a GridSearch and is called when parameters
#are estimated using optimization (either nloptr or genoud).
ExtractUnambiguousOptimizationParameters<-function(overall.results=NULL,results.list=NULL,migrationArray,
	rm.n0=TRUE,longNames=FALSE,sortParameters=TRUE,sortModelsAIC=TRUE,nonparmCols=2,optimization="grid"){

	#Add parameters
	all.parameters <- unique(sort(sapply(migrationArray,
		MsIndividualParametersConversionToUnambiguous,unambiguous.only = TRUE)))
	results.all <- c()
	row.count = 0

	#Get models
	whichRep<-1
	for (models in overall.results$models){ #Get current models
		param.mappings <- MsIndividualParametersConversionToUnambiguous(migrationArray[[whichRep]],
			longNames=longNames)
		param.mappings.temp<-param.mappings

		#Get nLoptr parameters
		param.unique<-unique(param.mappings.temp[,2])[which(unique(param.mappings.temp[,2]) != "0")]
		#If there is only 1 n0multiplier, toss it, as nLoptr doesn't return this value
		if(length(grep("n0multiplier",param.unique,value=TRUE)) == 1){
			oneN0multi<-TRUE
			param.unique<-param.unique[grep("n0multiplier",param.unique,invert=TRUE)]
			param.mappings.temp[,2][which(param.mappings.temp[,2] == "n0multiplier_1")]<-"1"
		}
		if(optimization == "nloptr"){
			nLoptr.parameters<-exp(results.list[[whichRep]]$solution)
			nLoptr.AIC<-results.list[[whichRep]]$objective
		}
		if(optimization == "genoud"){
			nLoptr.parameters<-exp(results.list[[whichRep]]$par)
			nLoptr.AIC<-results.list[[whichRep]]$value
		}
		names(nLoptr.parameters)<-param.unique

		##DoSingleRowExtraction Stuff
		result.vector <- rep(NA, 2+length(all.parameters))
		names(result.vector) <- c("models", "AIC", all.parameters)
		for(parameter.index in 2:length(result.vector)) {
			matching.element <- which(param.mappings.temp[,1] == names(result.vector)[parameter.index])
			if(length(matching.element)==1) {
				if(param.mappings.temp[matching.element,2]=="0") {
					result.vector[parameter.index] <- 0
				}
				if(param.mappings.temp[matching.element,2]=="1") {
					result.vector[parameter.index] <- 1
				}
				if(param.mappings.temp[matching.element,2] != "0" && param.mappings.temp[matching.element,2] != "1") {
					result.vector[parameter.index] <- as.numeric(nLoptr.parameters[param.mappings.temp[matching.element,2]])
				}
			}
		}
		result.vector[1]<-models
		if(optimization == "nloptr"){
			result.vector[2]<-results.list[[whichRep]]$objective
		}
		if(optimization == "genoud"){
			result.vector[2]<-results.list[[whichRep]]$value
		}
		results.all<-rbind(results.all,result.vector)
		whichRep<-whichRep + 1
	}

	#Convert all zero values to NA
	#Note that if tau = NA, this means tau was set to zero in the analysis and thus should stay
	#that way for model averaging
	results.all<-data.frame(results.all, check.names=FALSE, row.names=c(1:nrow(results.all)))
	if(longNames==TRUE){
		results.all[, !grepl("collapse",colnames(results.all))][,-(0:nonparmCols)][results.all[,
			!grepl("collapse",colnames(results.all))][,-(0:nonparmCols)] == 0] <- NA
	}else{
		results.all[, !grepl("t",colnames(results.all))][,-(0:nonparmCols)][results.all[,
			!grepl("t",colnames(results.all))][,-(0:nonparmCols)] == 0] <- NA
	}
	#Remove undesired parameters
	if(rm.n0==TRUE){
		if(longNames==TRUE){
    		results.all<-results.all[,!grepl("n0multiplier", colnames(results.all))]
    	}else{
    		results.all<-results.all[,!grepl("n", colnames(results.all))]
    	}
	}

	#Remove parameters that only contain NA's
	results.all <- RemoveParameterNAs(results.all)

	#Sort parameter columns into desired order
	if(sortParameters == TRUE){
		results.all<-sortParameterTable(results.all,(nonparmCols + 1),longNames=longNames)
	}

    #Sort models
    if(sortModelsAIC == TRUE){
    	results.all<-results.all[order(results.all$AIC),]
    }
    return(results.all)
}

#This function concatenates together results from separate phrapl runs using the merge function.
#It accommodates the new way that GridSearch calculates model averaged parameter values for unambiguous
#parameters (implemented in December 2015). These parameters are calculated in a GridSearch using the
#function ExtractUnambiguousGridParameters for Grid runs or the ExtractUnambiguousNLoptrParameters function
#for optimization runs.
#Just point to a directory of one or more rda files that can each contain phrapl results for one or more models.
#With older phrapl run output files, one can get the new unambiguous parameters using the function
#ConcatenateResults_unambiguousParametersForOldRuns.
#Or, to get ambiguous parameters (the old way) from these older result files, use the function
#ConcatenateResults_ambiguousParameters.
ConcatenateResults<-function(rdaFilesPath="./",rdaFiles=NULL, migrationArray, rm.n0=TRUE,
	longNames=FALSE, addTime.elapsed=FALSE, addAICweights=TRUE, outFile=NULL, nonparmCols=5){

    #If a vector of rda file names is not provided, then read them in from the given path
	if(is.null(rdaFiles)){
		rdaFiles<-grep(".rda",list.files(rdaFilesPath),value=TRUE)
	}

    #Stuff from ConcatenateResults
    totalData<-data.frame()
	result <- c()
	time.elapsed <- c()
	row.counter<-1
	templateRep<-c()
	if(addTime.elapsed == TRUE){
		nonparmColsWithTime<-nonparmCols + 1
	}
	nonparmColsOriginal<-nonparmCols

	#Extract and combine results across models
	for(rep in 1:length(rdaFiles)){
		load(paste(rdaFilesPath,rdaFiles[rep],sep="")) #Read model objects
		
		#Combine results from the rda into dataframe
		#If a result exists...
		if(!is.null(result)){
			continue=TRUE
			templateRep<-append(templateRep,rep) #Save those rep numbers that yield a result
			current.results<-result$overall.results[order(result$overall.results$models),] #sort by model number (same order as parameters)

 		#If a result is null for this model...
 		}else{
 		warning(paste("Warning: File",rdaFiles[rep],"contains no results.", sep=" "))
 		continue<-FALSE
 		
#			#Fill in null results for the null model (i.e., NAs)
#			current.results<-data.frame(matrix(NA,nrow=1,ncol=nonparmColsOriginal))
#			colnames(current.results)<-c("models","AIC","lnL","params.K","params.vector")
#			#First fill in the model info
#			current.results$models<-model
#			current.results$params.K<-KAll(migrationArray[[model]])
#			if(length(grep("n0multiplier",MsIndividualParameters(migrationArray[[model]]))) == 1){ #Toss n0multi, if there's one
#				current.results$params.vector<-paste(MsIndividualParameters(migrationArray[[model]])[-grep("n0multiplier",
#					MsIndividualParameters(migrationArray[[model]]))],collapse=" ")
#			}else{
#				current.results$params.vector<-paste(MsIndividualParameters(migrationArray[[model]]),collapse=" ")
#			}
#			#Then get the parameters
#			param.mappings <- MsIndividualParametersConversionToUnambiguous(migrationArray[[model]])
#			parametersNullNames<-param.mappings[which(param.mappings[,2] != "0"),][,1]
#			parametersNull<-data.frame(matrix(NA,nrow=1,ncol=length(parametersNullNames)))
#			colnames(parametersNull)<-parametersNullNames
#			current.results<-cbind(current.results,parametersNull)
		}

		#Add time if desired
		if(addTime.elapsed==TRUE & continue == T){
			load(paste(rdaFilesPath,rdaFiles[rep],sep="")) #make sure current model is loaded again

			if(length(current.results[,(nonparmColsOriginal + 1):ncol(current.results)]) == 1){
				singleParmName<-names(current.results[ncol(current.results)])
				current.results<-cbind(current.results[,1:nonparmColsOriginal],elapsedHrs=time.elapsed[,2],
					current.results[,(nonparmColsOriginal + 1):ncol(current.results)])
				colnames(current.results)[ncol(current.results)]<-singleParmName
			}else{
				current.results<-cbind(current.results[,1:nonparmColsOriginal],elapsedHrs=time.elapsed[,2],
					current.results[,(nonparmColsOriginal + 1):ncol(current.results)])
				nonparmCols<-nonparmColsWithTime
			}
		}

		#Add current.results to totalData
		#Note that with "merge", exact duplicate rows will not be merged. Thus, I have added a counting row to make each row
		#unique such that they can be merged. The row is tossed eventually.
		if(continue == T){
			row.col<-data.frame(matrix(row.counter:(row.counter + (nrow(current.results) - 1)),nrow=nrow(current.results),ncol=1))
			colnames(row.col)<-"row.col"
			current.results<-cbind(row.col,current.results)

			if(rep == 1){
				totalData<-rbind(totalData,current.results)
			}else{
				totalData<-merge(totalData,current.results,all=TRUE)
			}
			row.counter<-nrow(totalData) + 1
		}
	}
	
	#Remove undesired parameters
	if(rm.n0==TRUE){
		if(longNames==TRUE){
			if(length(grep("0multiplier", colnames(totalData))) <= 1){
    			totalData<-totalData[,!grepl("n0multiplier", colnames(totalData))]
    		}
    	}else{
    		if(length(grep("n", colnames(totalData))) <= 1){
    			totalData<-totalData[,!grepl("n", colnames(totalData))]
    		}
    	}
	}

	#Sort parameter columns into desired order
	totalData<-sortParameterTable(totalData,(nonparmCols + 2),longNames=longNames)

	#Make sure totalData is all sorted by model, then by AIC
	totalData<-totalData[order(totalData$models,totalData$AIC),]

	#Remove row counter column
	totalData<-totalData[,-1]
	row.names(totalData)<-1:nrow(totalData)

	#Add model ranks, dAIC, and wAIC
	if(addAICweights==TRUE){
		totalDataRank<-totalData[which(!is.na(totalData$AIC)),] #Break up dataset by AICs with and without NA
		totalDataNA<-totalData[which(is.na(totalData$AIC)),]
		dAICRankWeights<-GetAICweights(totalDataRank)
		totalDataRank<-cbind(totalDataRank[,1:3],dAICRankWeights,totalDataRank[,4:ncol(totalDataRank)])
		if(nrow(totalDataNA) > 0){
			totalDataNA<-cbind(totalDataNA[,1:3],dAICRankWeights[1:nrow(totalDataNA),],totalDataNA[,4:ncol(totalDataNA)])
			totalDataNA[,4:6]<-NA
			totalData<-rbind(totalDataRank,totalDataNA)
		}else{
			totalData<-totalDataRank
		}
		totalData<-totalData[order(totalData$dAIC,totalData$models),]
	}

	#Print table
	if(!is.null(outFile)){
		write.table(totalData,outFile,quote=FALSE,sep="\t",row.names=FALSE)
	}

	return(totalData)

}

#This function concatenates together results from separate phrapl runs.
#This differs from ConcatenateResults in that it can be used on phrapl results produced prior to adding
#unambiguous parameter estimation (but only for a grid, not if optimization was used).
#Thus, for old phrapl result files (prior to ~December 2015), this can model
#average parameter values directly from the AIC.Grid to produce unambiguous parameter estimates.
ConcatenateResults_unambiguousParametersForOldRuns<-function(rdaFilesPath="rdaFiles/",rdaFiles=NULL, migrationArray,
	rm.n0=TRUE, longNames=FALSE, nonparmCols=2, addTime.elapsed=FALSE, addAICweights=TRUE, outFile=NULL){

    #If a vector of rda file names is not provided, then read them in from the given path
	if(is.null(rdaFiles)){
		rdaFiles<-grep(".rda",list.files(rdaFilesPath),value=TRUE)
	}

    #Stuff from ConcatenateResults
    totalData<-data.frame()
	result <- c()
	time.elapsed <- c()
	row.counter<-1
	for (rep in 1:length(rdaFiles)){
		load(paste(rdaFilesPath,rdaFiles[rep],sep="")) #Read model objects

		#Combine results from the rda into dataframe
		current.results<-result$overall.results[order(result$overall.results$models),] #sort by model number (same order as parameters)
		if(addTime.elapsed==TRUE){
			current.results<-cbind(current.results,elapsedHrs=time.elapsed[,2])
		}

		#Get unambiguous, model-averaged parameters
		modelAveragedResults<-ExtractUnambiguousGridParameters(overall.results=result$overall.results,
			gridList=result$AIC.Grid,migrationArray=migrationArray[result$overall.results$models],
			nonparmCols=nonparmCols,longNames=longNames,rm.n0=rm.n0,sortParameters=FALSE,sortModelsAIC=FALSE)

    	#Add parameters estimates to current.results
    	current.results<-cbind(current.results,modelAveragedResults[-(0:nonparmCols)])

		#Add current.results to totalData
    	#Note that with "merge", exact duplicate rows will not be merged. Thus, I have added a counting row to make each row
    	#unique such that they can be merged. The row is tossed eventually.
		row.col<-data.frame(matrix(row.counter:(row.counter + (nrow(current.results) - 1)),nrow=nrow(current.results),ncol=1))
    	colnames(row.col)<-"row.col"
    	current.results<-cbind(row.col,current.results)

    	if(rep == 1){
			totalData<-rbind(totalData,current.results)
    	}else{
    		totalData<-merge(totalData,current.results,all=TRUE)
    	}
    	row.counter<-nrow(totalData) + 1
	}

	#Sort parameter columns into desired order
	beginParm<-ncol(result$overall.results)
	if(addTime.elapsed == TRUE){
		beginParm<-beginParm + 1
	}
    totalData<-sortParameterTable(totalData,(beginParm + 1),longNames=longNames)

	#Make sure totalData is all sorted by model
	totalData<-totalData[order(totalData$models),]

	#Remove row counter column
	totalData<-totalData[,-1]
	row.names(totalData)<-1:nrow(totalData)

	#Add model ranks, dAIC, and wAIC
	if(addAICweights==TRUE){
		dAICRankWeights<-GetAICweights(totalData)
		totalData<-cbind(totalData[,1:3],dAICRankWeights,totalData[,4:ncol(totalData)])
		totalData<-totalData[order(totalData$dAIC,totalData$models),]
	}

	#Print table
	if(!is.null(outFile)){
		write.table(totalData,outFile,quote=FALSE,sep="\t",row.names=FALSE)
	}

    return(totalData)
}

#This concatenates phrapl results, and also adds to these dAIC, wAIC, and model ranks for a set of models.
#For phrapl results calculated prior to the switch-over to unambiguous parameter names (~December 2015),
#this function summarizes parameters the "old way", using ambiguous parameter names. To get summaries using
#unambiguous parameter names from old results files, use the function ConcatenateResults_unambiguousParametersForOldRuns
ConcatenateResults_ambiguousParameters<-function(rdaFilesPath="./",rdaFiles=NULL,outFile=NULL,addAICweights=TRUE,
		rmNaParameters=TRUE,addTime.elapsed=FALSE,optimization=FALSE){

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
		if(optimization==TRUE){
			overall.results<-result[[3]][order(result[[3]][,1]),] #sort by model number (same order as parameters)
			if(addTime.elapsed==TRUE){
				current.results<-cbind(overall.results,elapsedHrs=time.elapsed[,2],result[[4]],result[[5]])
			}else{
				current.results<-cbind(overall.results,result[[4]],result[[5]])
			}
		}else{
			overall.results<-result[[2]][order(result[[2]][,1]),] #sort by model number (same order as parameters)
			if(addTime.elapsed==TRUE){
				current.results<-cbind(overall.results,elapsedHrs=time.elapsed[,2],result[[3]],result[[4]])
			}else{
				current.results<-cbind(overall.results,result[[3]],result[[4]])
			}
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

#Function for sorting parameter table
#Use in conjunction with ExtractUnambiguousGridParameters and ConcatenateResults
sortParameterTable<-function(overall.results,parmStartCol,longNames=FALSE){
    parmsOnly<-data.frame(matrix(as.matrix(overall.results[,parmStartCol:ncol(overall.results)]),nrow=nrow(overall.results),
    	ncol=length(overall.results[,parmStartCol:ncol(overall.results)])))
    colnames(parmsOnly)<-colnames(overall.results)[parmStartCol:ncol(overall.results)]
    if(longNames==TRUE){
	    overall.results<-cbind(overall.results[,((parmStartCol - parmStartCol) + 1):(parmStartCol - 1)],
    		parmsOnly[, grep("collapse",colnames(parmsOnly))][, order(colnames(parmsOnly[,
    			grep("collapse",colnames(parmsOnly))]))],
    		parmsOnly[, grep("n0multiplier",colnames(parmsOnly))][, order(colnames(parmsOnly[,
    			grep("n0multiplier",colnames(parmsOnly))]))],
    		parmsOnly[, grep("growth",colnames(parmsOnly))][, order(colnames(parmsOnly[,
    			grep("growth",colnames(parmsOnly))]))],
    		parmsOnly[, grep("migration",colnames(parmsOnly))][, order(colnames(parmsOnly[,
    			grep("migration",colnames(parmsOnly))]))])
    }else{
    	#If there is only a single parameter, need to coerce back into a data.frame
    	if(length(parmsOnly[1, grep("t",colnames(parmsOnly))]) == 1){
    		parms.t<-data.frame(matrix(parmsOnly[, grep("t",colnames(parmsOnly))]))
    		colnames(parms.t)<-colnames(parmsOnly)[grep("t",colnames(parmsOnly))]
    	}else{
    		parms.t<-parmsOnly[, grep("t",colnames(parmsOnly))][, sort(colnames(parmsOnly)[grep("t",colnames(parmsOnly))])]
    	}
    	if(length(parmsOnly[1, grep("n",colnames(parmsOnly))]) == 1){
    		parms.n<-data.frame(matrix(parmsOnly[, grep("n",colnames(parmsOnly))]))
    		colnames(parms.n)<-colnames(parmsOnly)[grep("n",colnames(parmsOnly))]
    	}else{
    		parms.n<-parmsOnly[, grep("n",colnames(parmsOnly))][, sort(colnames(parmsOnly)[grep("n",colnames(parmsOnly))])]
    	}
    	if(length(parmsOnly[1, grep("g",colnames(parmsOnly))]) == 1){
    		parms.g<-data.frame(matrix(parmsOnly[, grep("g",colnames(parmsOnly))]))
    		colnames(parms.g)<-colnames(parmsOnly)[grep("g",colnames(parmsOnly))]
    	}else{
    		parms.g<-parmsOnly[, grep("g",colnames(parmsOnly))][, sort(colnames(parmsOnly)[grep("g",colnames(parmsOnly))])]
    	}
    	if(length(parmsOnly[1, grep("m",colnames(parmsOnly))]) == 1){
    		parms.m<-data.frame(matrix(parmsOnly[, grep("m",colnames(parmsOnly))]))
    		colnames(parms.m)<-colnames(parmsOnly)[grep("m",colnames(parmsOnly))]
    	}else{
    		parms.m<-parmsOnly[, grep("m",colnames(parmsOnly))][, sort(colnames(parmsOnly)[grep("m",colnames(parmsOnly))])]
    	}

    	overall.results<-cbind(overall.results[,((parmStartCol - parmStartCol) + 1):(parmStartCol - 1)],
			parms.t,parms.n,parms.g,parms.m)
    }
    return(overall.results)
}

#Make a list of parameters for a model, renamed such that all parameters are unambiguous
#Used in conjunction with ExtractUnambiguousGridParameters
MsIndividualParametersConversionToUnambiguous<-function(migrationIndividual, unambiguous.only = FALSE,
	longNames=FALSE){
    collapseMatrix <- migrationIndividual$collapseMatrix
    complete <- migrationIndividual$complete
    n0multiplierMap <- migrationIndividual$n0multiplierMap
    growthMap <- migrationIndividual$growthMap
    migrationArray <- migrationIndividual$migrationArray
    parameterList <- c()
    unambiguousParameterList <- c()

	#First, assign population names for each generation
    populationNamesAssignments <- matrix(nrow = dim(collapseMatrix)[1],
        ncol = 1)
    populationNamesAssignments[, 1] <- sequence(dim(collapseMatrix)[1])
    if (max(collapseMatrix, na.rm = TRUE) > 0) {
        for (i in sequence(dim(collapseMatrix)[2])) {
            populationNamesAssignments <- cbind(populationNamesAssignments,
                populationNamesAssignments[, dim(populationNamesAssignments)[2]])
            merged <- sort(which(collapseMatrix[, i] == 1))
            for (merged.index in sequence(length(merged))) {
                populationNamesAssignments[merged[merged.index],
                  i + 1] <- populationNamesAssignments[merged[1],
                  i]
            }
            populationNamesAssignments[, i + 1] <- RenumberWithoutGaps(populationNamesAssignments[,
                i + 1])
        }
    }

	#Then define ancestral population names such that component tip population names are incorporated
    populationNameMatrix <- matrix(nrow = dim(populationNamesAssignments)[1],
        ncol = dim(populationNamesAssignments)[2])
    populationNameMatrix[, 1] <- as.character(populationNamesAssignments[,
        1])
    if (dim(populationNamesAssignments)[2] >= 2) {
        for (i in 2:dim(populationNamesAssignments)[2]) {
            for (j in sequence(dim(populationNamesAssignments)[1])) {
                populationNameMatrix[j, i] <- paste(sort(populationNameMatrix[which(populationNamesAssignments[,
                  i] == populationNamesAssignments[j, i]), i -
                  1]), collapse = "-")
                populationNameMatrix[j, i] <- paste(sort(unique(strsplit(populationNameMatrix[j,
                  i], "-")[[1]])), collapse = "-")
            }
        }
    }

    #Now with these new names, name collapse parameters
    if (max(collapseMatrix, na.rm = TRUE) > 0) {
        for (i in GetCollapseGenerationNamesWhenAddedEventsExist(collapseMatrix)) {
            parameterList <- paste("collapse_",GetCollapseGenerationNamesWhenAddedEventsExist(collapseMatrix),sep="")
            if(longNames == TRUE){
            	unambiguousParameterList <- append(unambiguousParameterList,
             	   paste("collapse_pop_", paste(populationNameMatrix[which(collapseMatrix[,
                	  i] >= 1), i], collapse = "_with_"), "_at_time_interval_",
                	  i, sep = ""))
            }else{
	           	unambiguousParameterList <- append(unambiguousParameterList,
             	   paste("t", i, "_", paste(populationNameMatrix[which(collapseMatrix[,
                	  i] >= 1), i], collapse = "."), sep = ""))
            }
        }
    }

    #Name n0 parameters
    for (row.index in sequence(dim(n0multiplierMap)[1])) {
        for (col.index in sequence(dim(n0multiplierMap)[2])) {
            focaln0value <- n0multiplierMap[row.index, col.index]
            if (!is.na(focaln0value)) {
                parameterList <- append(parameterList, paste("n0multiplier_",
                  focaln0value, sep = ""))
            	if(longNames == TRUE){
					unambiguousParameterList <- append(unambiguousParameterList,
					paste("n0multiplier_pop_", populationNameMatrix[row.index,
                    col.index], "_at_time_interval_", col.index,sep = ""))
            	}else{
                	unambiguousParameterList <- append(unambiguousParameterList,
						paste("n", col.index, "_", populationNameMatrix[row.index,
                   	 	col.index], sep = ""))
                }
            }
        }
    }

    #Name growth parameters
    for (row.index in sequence(dim(growthMap)[1])) {
        for (col.index in sequence(dim(growthMap)[2])) {
            focal.growth <- growthMap[row.index, col.index]
            if (!is.na(focal.growth)) {
            	if (focal.growth == 0){
            		parameterList <- append(parameterList, 0)
            	}else{
                	parameterList <- append(parameterList, paste("growth_",
                  		focal.growth, sep = ""))
                }
            	if(longNames == TRUE){
					unambiguousParameterList <- append(unambiguousParameterList,
					paste("growth_pop_", populationNameMatrix[row.index,
                    col.index], "_at_time_interval_", col.index,sep = ""))
            	}else{
                	unambiguousParameterList <- append(unambiguousParameterList,
						paste("g", col.index, "_", populationNameMatrix[row.index,
                   	 	col.index], sep = ""))
                }
            }
        }
    }

    #Name migration parameters
    for (from.index in sequence(dim(migrationArray)[1])) {
        for (to.index in sequence(dim(migrationArray)[2])) {
            for (time.index in sequence(dim(migrationArray)[3])) {
                focal.migration <- migrationArray[from.index,
                  to.index, time.index]
                if (!is.na(focal.migration)) {
					if (focal.migration == 0) {
						parameterList <- append(parameterList, 0)
					}else{
						parameterList <- append(parameterList, paste("migration_",
						focal.migration, sep = ""))
 					}
					if(longNames == TRUE){
						unambiguousParameterList <- append(unambiguousParameterList,
                    		paste("migration_from_", populationNameMatrix[from.index,
							time.index], "_to_", populationNameMatrix[to.index,
                      		time.index], "_at_time_interval_", time.index, sep = ""))
                    }else{
                    	unambiguousParameterList <- append(unambiguousParameterList,
                    		paste("m", time.index, "_", populationNameMatrix[from.index,
							time.index], ".", populationNameMatrix[to.index,
                      		time.index], sep = ""))
                    }
                }
            }
        }
    }

    return.object <- unambiguousParameterList
    if (!unambiguous.only) {
        return.object <- cbind(unambiguous = unambiguousParameterList,
            default = parameterList)
    }
    return(return.object)
}

#Extract parameters from a PHRAPL result for a given migrationIndividual
#To be used in conjunction with ExtractUnambiguousGridParameters
DoSingleRowExtraction <- function(x, param.mappings, all.parameters, models) {
	result.vector <- rep(NA, 2+length(all.parameters))
	names(result.vector) <- c("models", "AIC", all.parameters)
	for(parameter.index in 2:length(result.vector)) {
		matching.element <- which(param.mappings[,1] == names(result.vector)[parameter.index])
		if(length(matching.element)==1) {
			if(param.mappings[matching.element,2]=="0") {
				result.vector[parameter.index] <- 0
			} else {
				result.vector[parameter.index] <- as.numeric(x[param.mappings[matching.element,2]])
			}
		}
	}
	result.vector[1]<-models
	result.vector[2]<-as.numeric(x["AIC"])
	return(result.vector)
}

#Utility function used in conjunction with MsIndividualParametersConversionToUnambiguous
#Handles bookkeeping (stitches together missing entries: c(4, 1, 3, 1) becomes c(3, 1, 2, 1))
RenumberWithoutGaps<-function(x){
    return(as.numeric(as.factor(x)))
}

#This gets the generation times at which the collapses occur such that they match the
#generation times of other parameters, even when added events have been inserted into
#the collapseMatrix
GetCollapseGenerationNamesWhenAddedEventsExist<-function(collapseMatrix){
	#First, get the relative position of each event
	positionNames<-c(1:dim(collapseMatrix)[2])
	collapseGenerationNames<-c()
	for(i in 1:dim(collapseMatrix)[2]){
		if(sum(!is.na(collapseMatrix[,i])) > 0){
			if(max(collapseMatrix[,i],na.rm=TRUE) > 0){
				collapseGenerationNames<-append(collapseGenerationNames,i)
			}
		}
	}
	return(collapseGenerationNames)
}

#Function for getting model averaged parameter estimate for a particular model (across a grid)
#Use in conjunction with ExtractUnambiguousGridParameters
ModelAverageEachModel<-function(totalData,parmStartCol){
	first<-TRUE
	for(i in sort(unique(totalData$models))){
		totalData_m<-totalData[which(totalData$models == i),]
		totalData_m<-totalData_m[which(!is.na(totalData_m$AIC)),]
		totalData_m<-cbind(totalData_m[,((parmStartCol - parmStartCol) + 1):(parmStartCol - 1)],
			GetAICweights(totalData_m),totalData_m[parmStartCol:ncol(totalData_m)])
		modelAverages_m<-CalculateModelAverages(totalData_m,parmStartCol=(parmStartCol + 3))
		modelAverages_m<-cbind(totalData_m[1,((parmStartCol - parmStartCol) + 1):(parmStartCol - 1)],
			modelAverages_m)
		if(first == TRUE){
			modelAverages<-modelAverages_m
		}else{
			modelAverages<-merge(modelAverages,modelAverages_m,all=TRUE)
		}
		first<-FALSE
	}
	return(modelAverages)
}

#Calculate model averaged parameter values across a set of models
#The model averages can be calculated using one of two different methods, controlled by averageAcrossAllModels.
#If TRUE, each parameter is model averaged across all models in the dataset, even if that model
#does not include the given parameter. In these cases where a parameter is absent (i.e., the value
#is NA) this method simply assumes that the value for that parameter is zero. If FALSE, each
#parameter is model averaged by only considering those models that include the relevant parameter.

#The former method may be most appropriate for migration, as models that exclude this parameter
#are effectively setting migration to zero. However, the latter method may be more appropriate
#for the collapse time parameter (t) given that when a model excluding a particular coalescent
#event recieves high support, this can be subject to different interpretations. For example, this
#could signal that these two populations are so similar that coalescence time between them is
#effectively zero. However, it could also signal that these two populations are very distinct,
#but that they coalesce with other populations prior to coalescing with each other.

#parmStartCol gives the column number in totalData in which the parameter values begin. If using default
#PHRAPL output, this should be column 9.
CalculateModelAverages<-function(totalData,averageAcrossAllModels=TRUE,parmStartCol=9,keep.na=FALSE){
    #Add model averaged parameter estimates to a dataframe
    totalData <- totalData[order(totalData$AIC, totalData$models),]#sort parameters by AIC to match totalData
    totalData <- totalData[, grep(".*_I", colnames(totalData), invert = TRUE)]#Remove parameter index columns
    totalData <- totalData[which(!is.na(totalData$AIC)),] #Toss columns for which AIC is an NA

    #Remove parameter columns that contain only NAs
    if(keep.na == FALSE){
	    totalData <- RemoveParameterNAs(totalData)
	}

   	#If there is only a single model in the input, just return the values
   	if(nrow(totalData) == 1){
   		modelAverages<-totalData[,(parmStartCol:ncol(totalData))]
   		return(modelAverages)
   	}

   	#If model averaging across all models...
   	if(averageAcrossAllModels == TRUE){

	   #Calculate model averaged parameter estimates (NA's are ignored if present in some models)
		modelAverages <- data.frame(na.pass(totalData[, parmStartCol:ncol(totalData)] *
			totalData$wAIC))
		colnames(modelAverages)<-colnames(totalData)[-c(0:(parmStartCol - 1))]

		#Get model averages
		if(keep.na == TRUE){
			#If keeping parameter columns that contain only NAs, convert these model averages to NA...
			na.rm.dataset <- RemoveParameterNAs(modelAverages)
			NAcols<-which(colnames(modelAverages) %in% colnames(na.rm.dataset) == FALSE)
			modelAverages<-apply(modelAverages, 2, sum, na.rm = TRUE)
			modelAverages[NAcols]<-NA
		}else{
			modelAverages<-apply(modelAverages, 2, sum, na.rm = TRUE)
		}

		#Convert modelAverages into a dataframe
		modelAveragesMat <- data.frame(matrix(modelAverages, nrow = 1,
			ncol = length(names(modelAverages))))
		colnames(modelAveragesMat) <- names(modelAverages)
		modelAverages <- modelAveragesMat

	#Else, model average across only those models that contain the parameter
	}else{

  		modelAveragesVec<-c()
		for(j in parmStartCol:ncol(totalData)){
			totalDataTemp<-totalData[which(!is.na(totalData[,j])),]
			totalDataTemp<-totalDataTemp[order(totalDataTemp$models),] #GetAICweights exports in order of model
			wAIC.df<-GetAICweights(totalDataTemp)
			modelAverages<-totalDataTemp[,j] * wAIC.df$wAIC
			modelAverages<-sum(modelAverages)
			modelAveragesVec<-append(modelAveragesVec,modelAverages)
		}

		#If keeping parameter columns that contain only NAs, convert these model averages to NA...
		if(keep.na == TRUE){
			modelAveragesTemp<-data.frame(na.pass(totalData[, parmStartCol:ncol(totalData)] *
				totalData$wAIC))
			na.rm.dataset <- RemoveParameterNAs(modelAveragesTemp)
			NAcols<-which(colnames(modelAveragesTemp) %in% colnames(na.rm.dataset) == FALSE)
			modelAveragesVec<-apply(modelAveragesTemp, 2, sum, na.rm = TRUE)
			modelAveragesVec[NAcols]<-NA
		}

		modelAverages<-data.frame(matrix(modelAveragesVec,nrow=1,ncol=length(modelAveragesVec)))
		colnames(modelAverages)<-colnames(totalData)[parmStartCol:ncol(totalData)]
  	}
    return(modelAverages)
}

#Calculate dAICs, model ranks, and model weights
GetAICweights<-function(totalData=totalData){
	#Calculate change in AIC values from best model
	if(nrow(totalData) > 1){
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
	}else{
		#If there's only one row, of course don't need to model average
		dAICRankWeights<-data.frame("rank"=1,"dAIC"=0,"wAIC"=1)
		return(dAICRankWeights)
	}
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

#This function takes output from a phrapl run performed using multiple loci and re-calculates
#AICs and parameter estimates for different subsets of loci. This can either be done for independent
#subsets of loci (cumulative=FALSE) or for accumulating subsets of loci (cumulative=TRUE, in which
#the second subset is added to the first subset, and then the third subset is added to these, etc).
#lociRange gives the number of loci in the original analysis and NintervalLoci gives the desired number
#of loci in each subset, which must be a multiple of the total number of loci. If the models in
#the original analysis are only a subset of models in the migrationArray, this modelRange must be
#specified as well. Output of this function is an rda file for each locus subset (numbered in order
#with a numerical suffix, i.e., _1, _2, _3, etc). Note that subsets are taken in order from the original
#output files.
GenerateSetLoci<-function(lociRange,NintervalLoci,RoutFilename,rdaFilename,migrationArray,
	modelRange=c(1:length(migrationArray)),subsamplesPerGene,popAssignments,collapseStarts,
	migrationStarts,n0multiplierStarts,setCollapseZero=NULL,cumulative=FALSE,nTrees,
	dAIC.cutoff=2,nEq=nEq,totalPopVector=totalPopVector){

	#Get matching vectors from Rout
    matchesTEMP<-system(paste("grep matches -A1 ",RoutFilename," | grep [0123456789]",sep=""),intern=TRUE)
    matches<-c()
    for(i in 1:length(matchesTEMP)){
    	matches<-append(matches,strsplit(matchesTEMP[i],"\"")[[1]][2])
    }

	#Load rda file to get grid
	load(rdaFilename)

	##Sort GridList (sorted by AIC in rda, but want to be sorted by grid order)
	for(i in 1:length(gridList)){
		gridList[[i]]<-gridList[[i]][order(as.numeric(row.names(gridList[[1]]))),]
	}

    ##Creates a vector that repeats grid length per model
    sep_intervalsGrid<-c()
    for(modelID in 1:length(modelRange)){
        sep_intervalsGrid<-append(sep_intervalsGrid,c(as.numeric(rep(modelID,nrow(gridList[[modelID]])))))
    }

    ##Splits the matching vector into a list such that each model is an item in the list
    tableVectorByGridLength<-split(matches,sep_intervalsGrid)


	##################Print RDAs for each subset of loci

	#For each subset of loci...
	treesVec<-sequence(max(lociRange) * subsamplesPerGene)
	counterBegin<-1
	counterEnd<-NintervalLoci * subsamplesPerGene
	gridListOriginal<-gridList
	for(locusSubset in 1:(max(lociRange)/NintervalLoci)){ #Get positions of each locus subset in the matching vector
		gridList<-gridListOriginal
		currentLociRange<-treesVec[counterBegin:counterEnd]
		if(cumulative==FALSE){
			counterBegin<-counterBegin + (NintervalLoci * subsamplesPerGene)
		}
		counterEnd<-counterEnd + (NintervalLoci * subsamplesPerGene)

		#For each model...
		for(model in 1:length(migrationArray)){

			#For each parameter combination in the grid
            for(gridpoints in 1:length(tableVectorByGridLength[[model]])){
                totalVector<-as.numeric(strsplit(tableVectorByGridLength[[model]][gridpoints]," ")[[1]])
				totalVectorSubsampled<-totalVector[currentLociRange]

                lnLValue<-ConvertOutputVectorToLikelihood.1sub(outputVector=totalVectorSubsampled,
                	popAssignments=popAssignments,nTrees=nTrees,subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,
                	summaryFn="mean",nEq=nEq)

                #Replace the old AIC value in the grid with the new value (note: subtract from K for n0multiplier and any
                #collapses set to be zero.
                gridList[[model]][gridpoints,1]<-2*(-lnLValue + (KAll(migrationArray[[model]]) - length(setCollapseZero) - 1))
            }
        }

		#Re-sort gridList by AIC
		for(i in 1:length(gridList)){
			gridList[[i]]<-gridList[[i]][order(gridList[[i]]$AIC),]
		}

		#Now recalculate overall results object and parameter values object
        overall.results<-ExtractGridAICs(result=gridList,migrationArray=migrationArray,
        	modelRange=modelRange,setCollapseZero=setCollapseZero)
        parameters<-ExtractGridParameters(migrationArray=migrationArray,result=gridList,
			popVector=popAssignments[[1]],dAIC.cutoff=dAIC.cutoff)
		result<-list("AIC.Grid"=gridList,"overall.results"=overall.results,
			"parameters"=parameters[[1]],"parameterIndexes"=parameters[[2]])

		#Save all this to new rda - you get a new rda for each subset of loci
		save(list=c("gridList","overall.results","parameters","result","nTrees"),
			file=paste(rdaFilename,"_subset",locusSubset,".rda",sep=""))
	}
}

#The purpose of this function is to filter out high migration rates from the
#parameter grid after a phrapl GridSearch has been completed. The function
#inputs a gridList object, filters out migration in accordance with a specified
#criterion, and outputs the new gridList. One can then re-calculate overall AIC
#and parameter values for this gridList using ExtractGridAICs and
#ExtractGridParameters functions. If criterion = FilterMigrationGreaterThanCollapse,
#this will filter out any migration rate value that is greater than
#FilterMigrationGreaterThanCollapseMultiplier * the relevant collapse time. Thus, if
#FilterMigrationGreaterThanCollapseMultiplier = 1, this will filter out any
#parameter combination where migration is greater than tau.
FilterGridForMigration<-function(gridList,migrationArray,popAssignments,criterion=NULL,
	FilterMigrationGreaterThanCollapseMultiplier=1){

	for(whichModel in 1:length(gridList)){ #for each model

		currentParamNames<-MsIndividualParameters(migrationArray[[whichModel]])

		#Get vector the gives the migration matrix at which each migration parameter first occurs
		numMatrices<-length(migrationArray[[whichModel]]$migrationArray) /
			(length(popAssignments) * length(popAssignments))
		whichMigMatrixVec<-c() #to store the matrix at which each migration parameter first occurs
		for(j in 1:length(grep("migration",currentParamNames))){ #for each migration parameter
			for(i in 1:numMatrices){ #for each matrix
				if(length(grep(j,migrationArray[[whichModel]]$migrationArray[,,i])) > 0){
					whichMigMatrixVec<-append(whichMigMatrixVec,i)
					break()
				}else{}
			}
		}

		#Now filter the grid according to desired criterion
		gridList[[whichModel]]<-FilterMigrationByCriterion(gridListCurrent=gridList[[whichModel]],
			migrationIndividual=migrationArray[[whichModel]],whichMigMatrixVec=whichMigMatrixVec,
			criterion=criterion,FilterMigrationGreaterThanCollapseMultiplier=
			FilterMigrationGreaterThanCollapseMultiplier)
	}
	return(gridList)
}

#This function is called by FilterGridForMigration when filtering a gridList according to a specified criterion
FilterMigrationByCriterion<-function(gridListCurrent,migrationIndividual,whichMigMatrixVec,
	criterion=NULL,FilterMigrationGreaterThanCollapseMultiplier=1){
	paramNames<-MsIndividualParameters(migrationIndividual)
	for(i in 1:length(grep("migration",paramNames))){ #For each migration parameter
		currentCollapseColNum<-grep(paste("collapse_",whichMigMatrixVec[i],sep=""),colnames(gridListCurrent))
		currentMigrationColNum<-grep(paste("migration_",i,sep=""),colnames(gridListCurrent))

		#Criterion1: migration can't be larger than collapse
		if(criterion=="FilterMigrationGreaterThanCollapse"){
			gridListCurrent<-gridListCurrent[which((gridListCurrent[,currentCollapseColNum] *
				FilterMigrationGreaterThanCollapseMultiplier) >=
			gridListCurrent[currentMigrationColNum]),]
		}
	}
	return(gridListCurrent)
}

#This function calculates the genealogical divergence index using a given tau and migration rate
#from ms. To obtain this index, trees are iteratively simulated in ms with two species
#and three individuals under the given tau and migration rate (which can be asymmetrical).
#Then, the proportion of times that the given parameter values yield the true
#topology is calculated. This index gives something similar to the genealogical sorting index,
#except for a given set of migration and tau values instead of for groups on a tree (and for two
#focal taxa rather than for one focal taxon in relation to the rest of the tree): to what
#degree does the resulting location of alleles align with underlying species boundaries?
#We've scaled the index to be between 0 and 1, such that 0 = panmixia and 1 = strong support
#for divergence between whatever lineages are defined by the given tau and migration rate.
#Except when gdi is close to 1, there is a fairly high variance for the index, such that
#it should be calculated using a large number of replicates (10,000 at a minimum).
#It is possible to calculate a confidence interval. We have incorporated the binom.confint function
#from the binom package to do this, which can calculate a binomial confidence interval using
#eight different methods ("exact","asymptotic","agresti-coull","wilson","prop.test","bayes",
#"logit","cloglog","probit",or "profile"). Any of these can be specified using ciMethod, or "all"
#will result in all eight being calculated. If ciMethod=NULL, only the gdi will be outputted.
CalculateGdi<-function(tau,migration.in,migration.out=NULL,nreps=10000,ciMethod="exact",
	msPath=system.file("msdir","ms",package="P2C2M")) {
	if(is.null(migration.out)) {
		migration.out<-migration.in
	}
	divergence<-c()
	output<-system(paste(msPath," 3 ", nreps, " -T -I 2 2 1 -ej ",tau," 2 1 -m 1 2 ",migration.out," -m 2 1 ",migration.in,
		" | grep ';'",sep=""),intern=TRUE)
	intra.coalescence<-sapply(output,TestForWithinPopCoalescenceThreeSamplesOnly)
	#Return the proportion of three taxon trees where 1 and 2 coalesce before either coalesces with 3
	#Standardize so between 0 and 1
	if(is.null(ciMethod)){
		fraction<-sum(intra.coalescence)/length(intra.coalescence)
		fraction_standardized<-(fraction - 0.333) / (1 - 0.333)
		return(fraction_standardized)
	}else{
		ci<-binom.confint(x=sum(intra.coalescence),n=length(intra.coalescence),conf.level=0.95,
			methods=ciMethod)
		ci$mean<-(ci$mean - 0.333) / (1 - 0.333)
		ci$lower<-(ci$lower - 0.333) / (1 - 0.333)
		ci$upper<-(ci$upper - 0.333) / (1 - 0.333)
		return(ci)
	}
}

#Called within GetIntraspeciesCoalescentFraction
ConvertTreeString <- function(x) {
	return(read.tree(text=x))
}

#Called within GetIntraspeciesCoalescentFraction
TestForWithinPopCoalescenceThreeSamplesOnly <- function(x) {
	return(2==sum(grepl("1_2",GetClades(ConvertTreeString(x)))))
}

#The raw code as well as a binary executable for the program ms (Hudson 2002) comes pre-packaged with P2C2PM,
#a PHRAPL dependency. The purpose of \code{CheckMS} is to check whether the compiled binary file is compatible
#with one's operating system. If running this function returns a message that ms is working, then there is no
#need to compile ms prior to running PHRAPL. If the function returns an error or a message that there is
#something wrong, then you will need to compile ms using the \code{CompileMS} function. Once
#you have compiled ms, this check can again be run to ensure that the compilation was sucessful.
CheckMS<-function(){
	msPath<-system.file("msdir","ms",package="P2C2M")
	treeCheck<-suppressWarnings(system(paste(msPath," 1 1 -t 1",sep=""),intern=TRUE))
	if(length(grep(paste(msPath," 1 1 -t 1",sep=""),treeCheck[1])) > 0){
		return("ms is working")
	}else{
		return("Something is wrong. Run CompileMS")
	}
}

#If for whatever reason, the pre-packaged ms binary file doesn't work, you can compile the program
#on the fly using this function. ms can use one of two different random number generators.
#If compilation using the default \code{rand = 1} does not work, try \code{rand = 2}.
CompileMS<-function(rand=1){
	currentdir<-getwd()
	setwd(system.file("msdir",package="P2C2M"))
	if(rand==1){
		system("gcc -o ms ms.c streec.c rand1.c -lm")
	}else{
		system("gcc -o ms ms.c streec.c rand2.c -lm")
	}
	setwd(currentdir)
}

ReadStructureData <- function(file, pairs=2) {
	data <- read.delim(file=file, header=FALSE)
	sample.names <- as.character(data[,1])
	populations <- as.character(data[,2])
	data <- data[,c(-1, -2)]
	data[data==-9] <- NA
	rownames(data) <- sample.names
	if(pairs==2) {
		data.pruned <- matrix(nrow=length(sample.names), ncol=floor(dim(snps)[2]/2))
		for (i in sequence(dim(snps)[1])) {
			for (j in sequence(floor(dim(snps)[2]/2))) {
				data.pruned[i,j] <- data[i,sample(c(2*j, 2*j-1),1)]
			}
		}
		rownames(data.pruned) <- rownames(data)
		data <- data.pruned
	}
	if(pairs==1) {
		data.pruned <- matrix(nrow=length(sample.names)/2, ncol=dim(snps)[2])
		for (i in sequence(floor(dim(snps)[1]/2))) {
			for (j in sequence(dim(snps)[2])) {
				data.pruned[i,j] <- data[sample(c(2*i, 2*i-1),1),j]
			}
		}
		rownames(data.pruned) <- rownames(data)[seq(from=1, to=dim(data)[1]-1, by=2)]
		data <- data.pruned
	}
	return(list(snps=data, sample.names=sample.names, populations=populations))
}

ConvertStructureDataToTrees <- function(snps) {
	tree.vector <- c()
	for (i in sequence(dim(snps)[2])) {
		relevant.vector <- snps[,i]
		names(relevant.vector) <- rownames(snps)
		relevant.vector <- relevant.vector[!is.na(relevant.vector)]
		if (length(table(relevant.vector)) > 1 && min(table(relevant.vector))>1) { #need variation for a SNP
			left.side <- paste(names(relevant.vector)[which(relevant.vector==1)], collapse=",")
			right.side <- paste(names(relevant.vector)[which(relevant.vector==2)], collapse=",")
			tree.final <- ape::read.tree(text=paste('((', left.side, '),(', right.side, '));', sep=""))
			if(length(tree.vector)==0) {
				tree.vector <- c(tree.final)
			} else {
				tree.vector <- c(tree.vector, tree.final)
			}
		}
	}
	return(tree.vector)
}
