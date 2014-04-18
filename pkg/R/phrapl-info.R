#functions for getting or settiinfo about runs
SetMaxK<-function(popVector,nTrees=1,taxaPerParam=5,force=FALSE) {
  if (exists("maxK") && !force) {
    print(paste("note: setMaxK called, but maxK already exists, so using the existing one of ",maxK,". To change this behavior, give setMaxK the option force=TRUE",sep=""))
    return(maxK) 
  }
  else {
  	maxK<-max(1,floor(sum(popVector)/(nTrees*taxaPerParam))) 
  	return(maxK)
  }
}

MsIndividualParameters<-function(migrationIndividual) {
	collapseMatrix<-migrationIndividual$collapseMatrix
	complete<-migrationIndividual$complete
	n0multiplierMap<-migrationIndividual$n0multiplierMap
	migrationArray<-migrationIndividual$migrationArray
	parameterList<-c()
	if (max(collapseMatrix,na.rm=TRUE)>0) {
		for (i in sequence(KCollapseMatrix(collapseMatrix))) {
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



KAll<-function(migrationIndividual) {
  return(length(MsIndividualParameters(migrationIndividual)) -1) #-1 since first n0 is not free
}

KPopInterval<-function(popInterval) {
#returns the number of free parameters needed for that interval object. For example, having everything collapse in one step requires one param (the tMRCA). Having one collapse then a second requires two. Having no collapse requires 0
	maxVector<-ColMax(popInterval$collapseMatrix)
	return(length(which(maxVector>0)))
}

KCollapseMatrix<-function(collapseMatrix) {
#returns the number of free parameters needed for that interval object. For example, having everything collapse in one step requires one param (the tMRCA). Having one collapse then a second requires two. Having no collapse requires 0
	maxVector<-ColMax(collapseMatrix)
	return(length(which(maxVector>0)))
}


KN0multiplierMap<-function(n0multiplierMap) {
	return(max(n0multiplierMap,na.rm=TRUE)-1) #-1 since first n0 is not free
}



CreateAssignment<-function(popVector,file="assign.txt") {
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

CreateAssignment.df<-function(popVector) {
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


ReturnModel<-function(p,migrationArrayMap) {
   prunedResults<-subset(migrationArrayMap, migrationArrayMap$collapseMatrix.number==p[1])
   prunedResults<-subset(prunedResults, prunedResults$n0multiplierMap.number==p[2])
   prunedResults<-subset(prunedResults, prunedResults$migrationArray.number==p[3])
   if(dim(prunedResults)[1]==1) {
      return(prunedResults$model) 
   }
   else  {
   	  warning(paste("All these models match", prunedResults$model, "taking first one"))
   	  print(paste("All these models match", prunedResults$model, "taking first one"))
   	  print("Perhaps you sent ReturnModel a model ID rather than a vector?")
      return(prunedResults$model[1]) 
   }
}

ProportionDifference <- function(a,b) {
	return(abs((a-b)/min(c(a,b))))
}

PassBounds <- function(parameterVector, parameterBounds) {
	n0multiplierParameters<-sort(parameterVector[grep("n0multiplier",names(parameterVector))])
	migrationParameters<-sort(parameterVector[grep("migration",names(parameterVector))])
	collapseParameters<-sort(parameterVector[grep("collapse",names(parameterVector))])
	if(min(migrationParameters) < parameterBounds$minMigrationRate) {
		return(FALSE)
	}
	if(min(collapseParameters) <  parameterBounds$minCollapseTime) {
		return(FALSE)
	}
	if(length(n0multiplierParameters)>1) { #have at least two
		for (i in sequence(length(n0multiplierParameters)-1)) {
			if (ProportionDifference(n0multiplierParameters[i], n0multiplierParameters[i+1]) < parameterBounds$minN0Ratio) {
				return(FALSE)
			}
		}
	}
	if(length(migrationParameters)>1) { #have at least two
		for (i in sequence(length(migrationParameters)-1)) {
			if (ProportionDifference(migrationParameters[i], migrationParameters[i+1]) < parameterBounds$minMigrationRatio) {
				return(FALSE)
			}
		}
	}
	if(length(collapseParameters)>1) { #have at least two
		for (i in sequence(length(collapseParameters)-1)) {
			if (ProportionDifference(collapseParameters[i], collapseParameters[i+1]) < parameterBounds$minCollapseRatio) {
				return(FALSE)
			}
		}
	}
	if(length(collapseParameters)>1) { #have at least two
		for(i in sequence(length(collapseParameters)-1)){
			if(collapseParameters[i] > collapseParameters[i+1]){
				return(FALSE)
			}
		}	
	}
	return(TRUE)
}

#Return AIC for a given model and tree
ReturnAIC<-function(par,migrationIndividual,nTrees=1,msLocation="/usr/local/bin/ms",observed="observed.txt",print.results=FALSE, print.ms.string=FALSE, debug=FALSE, badAIC=100000000000000, maxParameterValue=100, parameterBounds=list(minCollapseTime=0.1, minCollapseRatio=0, minN0Ratio=0.1, minMigrationRate=0.05, minMigrationRatio=0.1),subsamplesPerGene=1,totalPopVector,popAssignments) {
  parameterVector<-exp(par)
  #now have to stitch in n0 being 1, always, for the first population
  positionOfFirstN0 <- min(grep("n0multiplier", MsIndividualParameters(migrationIndividual)))
  parameterVectorFirstPart<-parameterVector[sequence(positionOfFirstN0-1)]
  parameterVectorSecondPart<-parameterVector[(1+length(parameterVectorFirstPart)):length(parameterVector)]
  if((1+length(parameterVectorFirstPart)) > length(parameterVector)) {
  	parameterVectorSecondPart<-c()
  }
  parameterVector<-c(parameterVectorFirstPart, 1, parameterVectorSecondPart)
  names(parameterVector)<-MsIndividualParameters(migrationIndividual)
  if(!PassBounds(parameterVector, parameterBounds)) {
  	return(badAIC)
  }
  if(max(parameterVector)>maxParameterValue) {
  	return(badAIC)
  }
  if(debug) {
    print(parameterVector) 
  }
  likelihoodVector<-MatchingTrees(migrationIndividual=migrationIndividual,parameterVector=parameterVector,nTrees=nTrees,subsamplesPerGene=subsamplesPerGene,popAssignments=popAssignments,msLocation=msLocation,observed=observed,print.ms.string=print.ms.string, debug=debug)
  lnLValue<-ConvertOutputVectorToLikelihood(outputVector=likelihoodVector,nTrees=nTrees,subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,popAssignments=popAssignments)
  AICValue<-2*(-lnL.mat[1] + KAll(migrationIndividual))
  colnames(AICValue)<-"AIC"
	if(print.results) {
  		resultsVector<-c(AICValue,lnLValue[1],parameterVector)
    	names(resultsVector)<-c("AIC","lnL",MsIndividualParameters(migrationIndividual))
	    print(resultsVector)
 	  	print(paste(likelihoodVector,collapse=" ",sep="")) #print matches for each observed tree

		#print total number of matches
# 	   matches<-sum(as.numeric(likelihoodVector))
# 	   names(matches)<-"matchSum"
# 	   print(matches)

		#print number of matches per locus (summarized across subsamples)
#	    vector2print<-as.numeric(likelihoodVector)
#  		matches<-0
#    	matchesVec<-array()
#		localVector<-rep(NA, subsamplesPerGene)
#		print("matches per locus")
#		baseIndex<-1
#		for (i in sequence(length(vector2print))) {
#			localVector[baseIndex]<-vector2print[i]
#			baseIndex <- baseIndex+1
#			if(i%%subsamplesPerGene == 0) {
#				if(summaryFn=="SumDivScaledNreps"){
#					matchesVec<-rbind(matchesVec,get(summaryFn)(localVector=localVector,totalPopVector=totalPopVector,subNum=4,subsamplesPerGene=subsamplesPerGene,whichSampSize=min))
#				}else{
#					matchesVec<-rbind(matchesVec,get(summaryFn)(localVector))
#				}
#				localVector<-rep(NA, subsamplesPerGene)
#				baseIndex<-1
#			}
#		}
 #   	print(matchesVec[-1])
	}
  	return(cbind(AICValue,lnLValue))
}

#This takes an outputted grid from initial.AIC search and produces ranges of values for each parameter that can be constructed into
#a new, finer-grained grid. This new grid is based on the range of values contained in the best supported combinations of parameter values
#Outputted is a list of vectors containing a range of parameter values.
CreateFineGrid<-function(gridList=NULL,gridSizeVector=c(6,6,6)){
	#Create ranges of parameter values from which can construct new finer-grained grid
	gridSorted<-gridList[order(gridList$AIC),] #Sort current grid by AIC
	parmNames<-colnames(gridSorted)[-1]
	parmRanges<-data.frame(matrix(nrow=2,ncol=length(parmNames)))
	parmRangesExtendedTEMP<-data.frame(matrix(nrow=2,ncol=length(parmNames)))
	parmRangesExtended<-data.frame(matrix(nrow=2,ncol=length(parmNames)))
	fineGrid<-list()
	colnames(parmRanges)<-parmNames
	colnames(parmRangesExtendedTEMP)<-parmNames
	colnames(parmRangesExtended)<-parmNames
	for(rep in 1:length(parmNames)){ #for each parameter
		#Define ranges of the top 5 models
		parmRanges[1,rep]<-min(gridSorted[1:5,(rep + 1)])
		parmRanges[2,rep]<-max(gridSorted[1:5,(rep + 1)])
		uniqueParmsVec<-unique(gridSorted[,(rep + 1)])[order(unique(gridSorted[,(rep + 1)]))] #vector of unique parms
									
		#If the minimum best parameter is not the smallest parameter in the grid, extend the range a bit (to 50% the distance
		#from the next parameter value). Ditto for the maximum best parameter.
		if(is.na(uniqueParmsVec[rev(order(uniqueParmsVec[which(uniqueParmsVec < parmRanges[1,rep])]))][1])){
			parmRangesExtendedTEMP[1,rep]<-parmRanges[1,rep]
		}else{
			parmRangesExtendedTEMP[1,rep]<-uniqueParmsVec[rev(order(uniqueParmsVec[which(uniqueParmsVec < parmRanges[1,rep])]))][1]
		}
		if(is.na(uniqueParmsVec[which(uniqueParmsVec > parmRanges[2,rep])][1])){
			parmRangesExtendedTEMP[2,rep]<-parmRanges[2,rep]
		}else{
			parmRangesExtendedTEMP[2,rep]<-uniqueParmsVec[which(uniqueParmsVec > parmRanges[2,rep])][1]	
		}
									
		#Using the above range extentions, define grid ranges for each parameter
		parmRangesExtended[1,rep]<-parmRanges[1,rep] - ((parmRanges[1,rep] - parmRangesExtendedTEMP[1,rep]) / 2) 
		parmRangesExtended[2,rep]<-parmRanges[2,rep] + ((parmRangesExtendedTEMP[2,rep] - parmRanges[2,rep]) / 2) 

		#Now create vectors of new grid values for each parameter
		#First, define the number of grid points based on the type of parameter
		if(length(grep("collapse",colnames(parmRangesExtended)[rep])) != 0){
			gridSize<-gridSizeVector[1]
		}
		if(length(grep("n0multiplier",colnames(parmRangesExtended)[rep])) != 0){
			gridSize<-gridSizeVector[2]
		}
		if(length(grep("migration",colnames(parmRangesExtended)[rep])) != 0){
			gridSize<-gridSizeVector[3]
		}
									
		#Calculate internal grid values	
		gridInterval<-(parmRangesExtended[2,rep] - parmRangesExtended[1,rep]) / (gridSizeVector[1] - 1)
		gridInternalVals<-array()
		currentValue<-parmRangesExtended[1,rep]
		for(rep1 in 1:(gridSizeVector[1] - 2)){
			gridInternalVals<-c(gridInternalVals,(currentValue + gridInterval))
			currentValue<-tail(gridInternalVals,n=1)
		}
		gridInternalVals<-gridInternalVals[!is.na(gridInternalVals)]
		#Concatenate min, max, and internal values and add to the grid list
		fineGrid[[length(fineGrid) + 1]]<-c(parmRangesExtended[1,rep],gridInternalVals,parmRangesExtended[2,rep])
		names(fineGrid)[length(fineGrid)]<-colnames(parmRangesExtended)[rep]
		startGrid<-list()
		startGrid[[1]]<-log(expand.grid(fineGrid))									
	}
	return(fineGrid)
}

#This takes vectors of parameter values and constructs a grid containing all possible combinations of parameter values, except
#non-sensical collapse time combinations (e.g., tau's at time 1 which are larger than tau's at time 2) are filtered out.
#This grid is outputted in the form of a list which can be used to obtain AIC values using SearchContinuousModelSpaceNLoptr. 
CreateStartGrid<-function(fineGrid){
	startGrid<-list()
	startGrid[[1]]<-expand.grid(fineGrid)
	howManyCollapses<-length(grep("collapse",names(startGrid[[1]])))
	if(howManyCollapses > 1){
		for(rep in 1:(howManyCollapses - 1)){
			startGrid[[1]]<-startGrid[[1]][which(startGrid[[1]][,rep] < startGrid[[1]][,(rep+1)]),]
		}
	}
	startGrid[[1]]<-log(startGrid[[1]]) #since we are optimizing in log space
	return(startGrid)
}

##Match simulated trees to observed trees and export vector of matches
MatchingTrees<-function(migrationIndividual,parameterVector,popAssignments,nTrees=1,msLocation="/usr/local/bin/ms",observed="observed.txt",
	subsamplesPerGene,debug=FALSE,print.ms.string=FALSE,ncores=1) {

	msCallInfo<-CreateMSstringSpecific(popAssignments[[1]],migrationIndividual,parameterVector,ceiling(nTrees/ncores))
	
  if(print.ms.string) {
    print(msCallInfo) 
  }
  if(debug) {
    print("parameterVector")
    print(parameterVector)
  }

	#Simulate trees in MS
	outputstringMS<-paste(msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts," | grep ';' > mstrees.txt", sep=" ") 
	system(outputstringMS,intern=TRUE)
	
	#Swap simulated individual labels with population labels
	phy<-read.tree("mstrees.txt")
	for(i in 1:length(phy)){
		phy[[i]]$tip.label<-ConvertAlleleNumbersToPopulationLetters(phy=phy[[i]],popVector=popAssignments[[1]])
	}	
	
	#Swap observed individual labels with population labels
	observed<-read.tree("observed.tre")
	treesPerLocus<-subsamplesPerGene * length(popAssignments) #number of trees per locus
	nLoci<-length(observed) / treesPerLocus #number of loci
	treeCounter<-1
	for(e in 1:nLoci){ #for each locus
		subsampSizeCounter<-1
		for(f in 1:(length(observed) / nLoci)){ #for each locus
			#Swap observed labels
			observed[[treeCounter]]$tip.label<-ConvertAlleleNumbersToPopulationLetters(phy=observed[[treeCounter]],popVector=popAssignments[[subsampSizeCounter]])
			treeCounter<-treeCounter + 1
			if(f%%subsamplesPerGene == 0){ #increase to next popAssignments size class
				subsampSizeCounter<-subsampSizeCounter + 1
			}
		}
	}

	#Calculate popAssignment indexes for the matching vector
	outputVector<-array(NA,length(observed))
	counter1<-1
	for(n in 1:nLoci){
		counter2<-1
		for(m in 1:(length(observed) / subsamplesPerGene)){
			outputVector[counter1]<-counter2
			counter1<-counter1+1	
			if(m%%subsamplesPerGene==0){
				counter2<-counter2+1
			}
		}
	}

	#Calculate number of matches for each observed tree
	currentPhy<-phy
	GetClades <- function(phy) {
		return(simplify2array(sapply(subtrees(phy), GetAndBindLabel)))
	}
	for(q in 1:length(popAssignments)){ #for each popAssignment size class
		for(r in which(outputVector==q)){ #for each subsample within a size class
			cladesAllPhy<-lapply(currentPhy, GetClades) #get clades for each simulated tree for this size class
			matches<-0
			cladesGene <- simplify2array(sapply(subtrees(observed[[r]]), GetAndBindLabel)) #get clades for all the observed trees
			for(s in 1:length(currentPhy)){ #for each simulated tree
				matches<-matches + GetScoreOfSingleTree(cladesMS=cladesAllPhy[[s]], phyMS=currentPhy[[s]],
					cladesGene=cladesGene, phyGene=observed[[r]])
			}
			outputVector[r]<-matches
		}
		if(q < length(popAssignments)){
			for(t in 1:length(currentPhy)){
				currentPhy[[t]]<-SubsampleMSTree(phy=currentPhy[[t]],popVectorOrig=popAssignments[[q]],popVectorFinal=popAssignments[[q+1]])
			}
		}
	}
	
	return(outputVector)

}

#This function calculates lnL based on number of matches for each subsample. Then the likelihood of the
#full dataset is estimated from the subsampled dataset based on a linear relationship between sample size
#and lnL.
ConvertOutputVectorToLikelihood<-function(outputVector,popAssignments,nTrees=1,subsamplesPerGene=1,totalPopVector){
	treesPerLocus<-subsamplesPerGene * length(popAssignments) #number of trees per locus
	nLoci<-length(outputVector) / treesPerLocus #number of loci
	
	#Adjust zero matches to equal probability of getting match
	beginLocus<-1
	whichLocus<-1
	outputVector<-as.numeric(outputVector)
	for(f in 1:nLoci){ #for each locus
		localVector<-outputVector[beginLocus:(whichLocus * treesPerLocus)]
		beginSizeClass<-1
		endSizeClass<-length(localVector) / length(popAssignments)
		for(g in 1:length(popAssignments)){ #for each subsample size class
			probOfMissing<-1/howmanytrees(sum(popAssignments[[g]]))
			localVector[beginSizeClass:endSizeClass][which(localVector[beginSizeClass:endSizeClass]==0)]<-probOfMissing
			beginSizeClass<-beginSizeClass + subsamplesPerGene
			endSizeClass<-endSizeClass + subsamplesPerGene
		}
		outputVector[beginLocus:(whichLocus * treesPerLocus)]<-localVector
		beginLocus<-beginLocus + treesPerLocus
		whichLocus<-whichLocus + 1
	}
	
	#Convert to log space
	outputVector<-outputVector/nTrees
	outputVector<-log(outputVector)
	
	#Calculate coefficients for ntaxa-lnL relationship to calculate lnL of full dataset
	beginLocus<-1
	whichLocus<-1
	slopes<-data.frame(matrix(NA,nrow=subsamplesPerGene,ncol=nLoci))
	intercepts<-data.frame(matrix(NA,nrow=subsamplesPerGene,ncol=nLoci))
	#For each locus	
	for(i in 1:nLoci){
		localVector<-outputVector[beginLocus:(whichLocus * treesPerLocus)]
		beginLocus<-beginLocus + treesPerLocus
		whichLocus<-whichLocus + 1
		currentSlopes<-data.frame()
		currentIntercepts<-data.frame()
		#For each subsample (but across different subsample size classes)
		for(j in 1:(length(localVector) / length(popAssignments))){
			#Get indexes to find a given subsampled tree of different sizes
			currentSubsample<-array(NA,length(popAssignments))
			currentSubsampleSize<-j
			for(k in 1:length(popAssignments)){
				currentSubsample[k]<-currentSubsampleSize
				currentSubsampleSize<-currentSubsampleSize + subsamplesPerGene
			}
			#Get the slopes and intercepts for each subsample within the current locus
			slopes[j,i]<-coef(lm(localVector[currentSubsample]~sapply(popAssignments,sum)))[2]
			intercepts[j,i]<-coef(lm(localVector[currentSubsample]~sapply(popAssignments,sum)))[1]
		}
		colnames(slopes)[i]<-paste("slopes.L",i,sep="")
		colnames(intercepts)[i]<-paste("intercepts.L",i,sep="")
	}
	#Get mean and 95%CI of coefficients
	meanSlopes<-sapply(slopes,mean,simplify = "array")
	meanIntercepts<-sapply(intercepts,mean,simplify = "array")	
	coeffsOrdered<-sapply(cbind(slopes,intercepts),sort)
	if(subsamplesPerGene>1){
		coeffs95Lower<-coeffsOrdered[(round(subsamplesPerGene*0.025) + 1),]
		coeffs95Upper<-coeffsOrdered[round(subsamplesPerGene*0.975),]
	}else{
		coeffs95Lower<-unname(coeffsOrdered)
		coeffs95Upper<-unname(coeffsOrdered)		
	}
	names(coeffs95Lower)<-paste("L95.",c(names(slopes),names(intercepts)),sep="")
	names(coeffs95Upper)<-paste("U95.",c(names(slopes),names(intercepts)),sep="")
	#Estimate lnL from full dataset
	lnL.byLocus<-unname(meanIntercepts + (sum(totalPopVector) * meanSlopes))
	names(lnL.byLocus)<-paste("lnL.L",c(1:length(lnL.byLocus)),sep="")
	lnL<-sum(lnL.byLocus)
	names(lnL)<-"lnL"
	#Prep to return
	lnL<-c(lnL,lnL.byLocus,meanSlopes,meanIntercepts,coeffs95Lower,coeffs95Upper)
	lnL.mat<-data.frame(matrix(nrow=1,ncol=length(lnL),lnL))
	colnames(lnL.mat)<-names(lnL)

	return(lnL.mat)
}

ConvertAlleleNumbersToPopulationLetters <- function(phy, popVector) {
	assignFrame <- CreateAssignment.df(popVector)
	phy$tip.label <- as.character(assignFrame$popLabel[match(phy$tip.label, assignFrame$indivTotal)])
}

#Drops tips from the simulated trees (to be used in accordance with popAssignments)
SubsampleMSTree <- function(phy, popVectorOrig, popVectorFinal) {
	taxaToDrop <- c()
	minSample <- 1
	maxSample <- 0
	uniquePops<-unique(phy$tip.label)
	for (population in 1:length(uniquePops)) {		
		maxSample <- minSample + popVectorOrig[population] - 1
		sampleSizeDiff <- popVectorOrig[population] - popVectorFinal[population]
		if (sampleSizeDiff > 0) {
			taxaToDrop <- append(taxaToDrop, 
				which(phy$tip.label==uniquePops[population])[c(1:sampleSizeDiff)])
		}
		minSample <- maxSample + 1
	}
	phy <- drop.tip(phy, taxaToDrop)

	return(phy)
}

#Gets a list of all clades in the tree. For the taxa descended from each clade, sorts alphabetically and then makes them a string
GetAndBindLabel <- function(phy) { 
	#note the sorting here
	return( paste( sort( phy$tip.label ), collapse="_" ) )
}

#Gets degree of polytomies in order to apply a correction when matching trees
GetOutDegreeOfPolytomies <- function(phy) {
	descendantCounts <- table(phy$edge[,1])
	descendantCounts <- unname(descendantCounts[which(descendantCounts>2)])
	return(descendantCounts)
}

#Match clades of two trees
#assumes:
#1. You have already run ConvertAlleleNumbersToPopulationLetters so these trees have letters
#2. You have already made them have the same size (do SubsampleMSTree if needed)
GetScoreOfSingleTree <- function(cladesMS, phyMS, cladesGene, phyGene) {
	#cladesMS <- simplify2array(sapply(subtrees(phyMS), GetAndBindLabel))
	#cladesGene <- simplify2array(sapply(subtrees(phyGene), GetAndBindLabel))
	numberCladesInMSOnly <- sum(!cladesMS%in%cladesGene)
	numberCladesInGeneOnly <- sum(!cladesGene%in%cladesMS)
	matchCount <- 0
	if(numberCladesInMSOnly == 0 && numberCladesInGeneOnly==0) {
		matchCount <- 1
	}
	if (numberCladesInMSOnly > 0 && numberCladesInGeneOnly == 0) {
		descendantCounts <- GetOutDegreeOfPolytomies(phyGene)
		correction <- 1
		for (i in sequence(length(descendantCounts))) {
			correction <- correction * howmanytrees(descendantCounts[i], rooted=TRUE, binary=TRUE, labeled=TRUE)
		}
		matchCount <- 1 / correction #idea here is that the gene tree could have been resolved at each polytomy multiple ways
		#only one of these ways would match the given phyMS tree. So we figure out the number of ways to resolve polytomy 1, multiply that by the number of ways to resolve polytomy 2, etc. A polytomy with three descendant edges has 3 ways to resolve it, one with 4 descendant edges has 3 * 5 = 15, etc.
	}
	return(matchCount)
}