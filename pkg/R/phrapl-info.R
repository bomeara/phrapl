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


FractionNonZeroMigration <- function(migrationIndividual) {
	return(sum(migrationIndividual$migrationArray[!is.na(migrationIndividual$migrationArray)]>0) / sum(!is.na(migrationIndividual$migrationArray)) )	
}

NumberPopulationsAtRoot <- function(migrationIndividual) {
	last.interval <- dim(migrationIndividual$collapseMatrix)[2]
	number.ones <- length(which(migrationIndividual$collapseMatrix[,last.interval]==1))	
	number.na <- length(which(is.na(migrationIndividual$collapseMatrix[,last.interval])))
	number.all <- dim(migrationIndividual$collapseMatrix)[1]
	number.alive <- number.all - number.na
	number.unmerged <- number.alive
	if (number.ones > 0) {
		number.unmerged <- number.alive - number.ones + 1
	}
	return(number.unmerged)
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

GuessAssignments <- function(observedTrees, popVector) {
	tips <- observedTrees[[1]]$tip.label
	# Levenshtein Distance
	distances  <- adist(tips)
	rownames(distances) <- tips
	hc <- hclust(as.dist(distances))
	groupings <- cutree(hc, k=length(popVector))
	if(sum(sort(table(groupings))==sort(popVector))!=length(popVector)) {
		warning("automatic grouping failed")
		return(NA)
	}
	alphabet<-strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ","")[[1]]
	assignments<-data.frame(indiv=names(groupings), popLabel=groupings, row.names=sequence(length(groupings)), stringsAsFactors=FALSE)
	groupings.table <- table(groupings)
	for(i in sequence(length(groupings.table))) {
		matching.index<-which(popVector==groupings.table[i])
		if(length(matching.index)>0) {
			matching.index<-matching.index[1]
		} else {
			warning("automatic grouping failed")
			return(NA)
		}
		assignments$popLabel[which(assignments$popLabel==i)]<-alphabet[matching.index]
		popVector[matching.index]<-0
		
	}
	return(assignments)
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
ReturnAIC<-function(par,migrationIndividual,nTrees=1,msPath="ms",comparePath=system.file("extdata", "comparecladespipe.pl", package="phrapl"),
		subsampleWeights.df=NULL,
		unresolvedTest=TRUE,print.results=FALSE, print.ms.string=FALSE,debug=FALSE,print.matches=FALSE,
		badAIC=100000000000000,ncores=1,maxParameterValue=100,parameterBounds=list(minCollapseTime=0.1,
		minCollapseRatio=0,minN0Ratio=0.1,minMigrationRate=0.05,minMigrationRatio=0.1),subsamplesPerGene=1,
		totalPopVector,popAssignments,summaryFn="mean",saveNoExtrap=FALSE, doSNPs=FALSE, nEq=100){
	parameterVector<-exp(par)
	#now have to stitch in n0 being 1, always, for the first population
	positionOfFirstN0 <- grep("n0multiplier", MsIndividualParameters(migrationIndividual))[1]
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
 	if(print.results) {
 	   print(parameterVector) 
 	}

 	#Do tree matching
 	likelihoodVector<-c()
	for(j in 1:length(popAssignments)){ #Do separately for each subsample size class
		currentPopAssign<-j
		likelihoodVectorCurrent<-PipeMS(migrationIndividual=migrationIndividual,parameterVector=parameterVector,
		nTrees=nTrees,subsamplesPerGene=subsamplesPerGene,popAssignments=popAssignments,msPath=msPath,comparePath=comparePath,
		ncores=ncores,currentPopAssign=currentPopAssign,print.ms.string=print.ms.string,debug=debug,unresolvedTest=unresolvedTest, doSNPs=doSNPs)

 	 	#Apply weights to matches
		if(!is.null(subsampleWeights.df)) {
			likelihoodVectorCurrent<-as.numeric(likelihoodVectorCurrent) * subsampleWeights.df[[j]][,1]
		} else {
			likelihoodVectorCurrent<-as.numeric(likelihoodVectorCurrent)
		}

	  	likelihoodVector<-append(likelihoodVector,likelihoodVectorCurrent)
  	}

  	#Convert matches to likelihoods
 	if(length(popAssignments) > 1){
  		lnLValue<-ConvertOutputVectorToLikelihood(outputVector=likelihoodVector,popAssignments=popAssignments,
  			nTrees=nTrees,subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,saveNoExtrap=saveNoExtrap,
  			summaryFn=summaryFn, nEq=nEq) 
  		AICValue<-2*(-lnLValue[1] + KAll(migrationIndividual))
  		colnames(AICValue)<-"AIC"
  		if(saveNoExtrap==TRUE){
  			AICValue.noExtrap<-2*(-lnLValue[2] + KAll(migrationIndividual))
  			colnames(AICValue.noExtrap)<-"AIC.lnL.noExtrap"
  		}		
  	}else{
  		lnLValue<-ConvertOutputVectorToLikelihood.1sub(outputVector=likelihoodVector,
  			popAssignments=popAssignments,nTrees=nTrees,subsamplesPerGene=subsamplesPerGene,
  			totalPopVector=totalPopVector,summaryFn=summaryFn, nEq=nEq)	
  		AICValue<-2*(-lnLValue[1] + KAll(migrationIndividual))		
	}
	if(print.results) {
		parameterVector<-as.data.frame(matrix(nrow=1,ncol=length(parameterVector),parameterVector))
		resultsVector<-cbind(AICValue,lnLValue[1],parameterVector)
    	names(resultsVector)<-c("AIC","lnL",MsIndividualParameters(migrationIndividual))
#	    print(resultsVector)
		print(resultsVector[1:2])
		if(print.matches){
 	  		cat("\nmatches\n")
 	  		print(paste(round(likelihoodVector,4),collapse=" ",sep="")) #print matches for each observed tree
 	  	}
		cat("End Run\n\n") #Separator of different simulation runs

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
#					matchesVec<-rbind(matchesVec,get(summaryFn)(localVector=localVector,totalPopVector=totalPopVector,subNum=4,subsamplesPerGene=subsamplesPerGene))
#				}else{
#					matchesVec<-rbind(matchesVec,get(summaryFn)(localVector))
#				}
#				localVector<-rep(NA, subsamplesPerGene)
#				baseIndex<-1
#			}
#		}
 #   	print(matchesVec[-1])
	}
	if(length(popAssignments) > 1 && saveNoExtrap==TRUE){
		return(cbind(AICValue,AICValue.noExtrap,lnLValue))
	}else{
  		return(cbind(AICValue,lnLValue))
  	}
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
PipeMS<-function(migrationIndividual,parameterVector,popAssignments,nTrees=1,msPath="ms",
		comparePath=system.file("extdata", "comparecladespipe.pl", package="phrapl"),unresolvedTest=TRUE,subsamplesPerGene,debug=FALSE,
		print.ms.string=FALSE,ncores=1,currentPopAssign=1, doSNPs=FALSE){
	stored.wd<-getwd()
	setwd(tempdir())
	observed<-paste(tempdir(),"/observed",currentPopAssign,".tre",sep="")
	assign<-paste(tempdir(),"/assign",currentPopAssign,".txt",sep="")
	msCallInfo<-CreateMSstringSpecific(popAssignments[[currentPopAssign]],migrationIndividual,
		parameterVector,ceiling(nTrees/ncores))
	
	systemMS<-function(stringname){
#		outputVectorMS<-suppressWarnings(system(command="(stringname) & sleep 2 ; kill $!",intern=TRUE))
#		http://stackoverflow.com/questions/687948/timeout-a-command-in-bash-without-unnecessary-delay
		outputVectorMS<-system(stringname,intern=TRUE)
		return(outputVectorMS)
	}
	systemPerl<-function(stringname){
		outputVectorPerl<-system(stringname,intern=TRUE)
		return(outputVectorPerl)
	}

	if(print.ms.string) {
    	print(msCallInfo) 
  	}

	unresolvedFlag<-"-u"
	if (unresolvedTest==FALSE) {
		unresolvedFlag<-""
	}

	snpFlag<-"-d" #for doSNPs
	if (doSNPs==FALSE) {
		snpFlag<-""
	}

	
	#Simulate trees in MS and do matching in perl
	outputstringMS<-paste(msPath,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts,
		" | grep ';' > mstrees.txt", sep=" ") 
	outputstringPerl<-paste("cat mstrees.txt | perl",comparePath,unresolvedFlag, snpFlag, paste("-a",assign,sep=""), 
		paste("-o",observed,sep=""),sep=" ")

	if(ncores==1){
		outputVectorMS<-systemMS(stringname=outputstringMS)
		outputVectorPerl<-systemPerl(stringname=outputstringPerl)
		outputVector<-paste(outputVectorMS,outputVectorPerl,sep=" ")
		setwd(stored.wd)
		return(outputVector)
	}else{
		stop("something is wrong here")
		wrapOutput<-function(x,outputstring) {
			as.numeric(system(outputstring,intern=TRUE))
		}
		outputVector<-apply(simplify2array(mclapply(sequence(ncores),wrapOutput,outputstring=paste(outputstringMS, outputstringPerl, sep=" ", collapse=" "),mc.cores=ncores)),1,sum)
		setwd(stored.wd)
		return(outputVector)
	}
}

#What happens with zero matching trees? Is the probability of the data exactly
#zero under that model? No -- there is a nonzero probability of nearly every
#gene tree topology, and letting it be zero is a strong statement that a model
#is impossible. It's like flipping a coin 5 times, finding zero heads, and
#saying that the best estimate for P(heads) is zero: it's probably wrong.
#Similarly, we have prior information: each possible gene topology, absent
#other information, has a 1/howmanytrees(nTips) chance of matching. Simply
#plugging this in when we have no matches isn't ideal: shouldn't our estimate
#of this change if we do 1e9 simulated trees than 1e3 and still find no match?
#We use a beta-binomial distribution to do this; mean is set to
#1/howmanytrees(nTips), and the equivalent sample size of the prior
#(how much weight we give it) matters but is up to the user
#By default, we assume that it's equal to 100 data points
#This is several orders of magnitude smaller than the usual number of nTrees
#so using beta will have little effect with matches, but it does let lack of
#matches result in a finite likelihood
AdjustUsingBeta <- function(numMatches, nTrees, nTips, nEq=100) {
	#nEq = a + b + 1
	#1/howmanytrees(nTips) = a / (a + b) #the prior mean
	#b = nEq - a - 1
	#howmanytrees(nTips) = (a + b) / a = (a + nEq - a - 1) / a = (nEq -1) / a
	a <- (nEq - 1) /  howmanytrees(nTips)
	b <- nEq - a - 1
	return( (numMatches + a) / (nTrees + a + b))	
}

#This function calculates lnL based on number of matches for each subsample. Then the likelihood of the
#full dataset is estimated from the subsampled dataset based on a linear relationship between sample size
#and lnL.
ConvertOutputVectorToLikelihood<-function(outputVector,popAssignments,nTrees=1,subsamplesPerGene=1,totalPopVector,
		saveNoExtrap=FALSE,summaryFn="mean", nEq=100){
	nLoci<-length(outputVector) / length(popAssignments) / subsamplesPerGene #number of loci
		
	start<-1
	end<-length(outputVector) / length(popAssignments)
	for(g in sequence(length(popAssignments))){ #for each subsample size class
		currentPopAssignments<-popAssignments[[g]]
		outputVector[start:end]<-as.numeric(outputVector[start:end])
		#probOfMissing<-1/((howmanytrees(sum(popAssignments[[g]]))) * 1000000)
		#outputVector[start:end][which(outputVector[start:end]==0)]<-probOfMissing
		outputVector[start:end]<-AdjustUsingBeta(numMatches=outputVector[start:end], nTrees=nTrees, nTips=sum(currentPopAssignments), nEq=nEq)

		start<-start + end
		end<-end + end
	}
	
	#Convert to log space
#	outputVector<-as.numeric(outputVector)/nTrees
	outputVector<-log(outputVector)

	#If the option to save non-extrapolated likelihood using the largest subsample size, do this
	if(saveNoExtrap==TRUE){
		lnL.noExtrap<-0
		localVector.noExtrap<-rep(NA,subsamplesPerGene)
		baseIndex<-1
		locusIndex<-1
		for (i in 1:(length(outputVector) / length(popAssignments))){
			localVector.noExtrap[baseIndex]<-outputVector[i]
			baseIndex <- baseIndex+1
			if(i%%subsamplesPerGene == 0) {
				if(summaryFn=="SumDivScaledNreps"){
					lnL.noExtrap<-lnL.noExtrap+get(summaryFn)(localVector=localVector.noExtrap,popAssignments,
						totalPopVector=totalPopVector[locusIndex],subsamplesPerGene=subsamplesPerGene)
				}else{
					lnL.noExtrap<-lnL.noExtrap+get(summaryFn)(localVector.noExtrap)
				}
				localVector.noExtrap<-rep(NA, subsamplesPerGene)
				baseIndex<-1
				locusIndex<-locusIndex + 1
			}
		}
	}
		
	#Index locations of loci within outputVector
	treesPerLocus<-subsamplesPerGene * length(popAssignments) #number of trees per locus
	outputIndex<-array(NA,length(outputVector))
	counter1<-1
	for(n in 1:length(popAssignments)){
		counter2<-1
		for(m in sequence(length(outputIndex) / length(popAssignments))){
			outputIndex[counter1]<-counter2
			counter1<-counter1+1	
			if(m%%subsamplesPerGene==0){
				counter2<-counter2+1
			}
		}
	}

	#Calculate coefficients for ntaxa-lnL relationship to calculate lnL of full dataset
	slopes<-data.frame(matrix(NA,nrow=subsamplesPerGene,ncol=nLoci))
	intercepts<-data.frame(matrix(NA,nrow=subsamplesPerGene,ncol=nLoci))	
	#For each locus	
	for(i in 1:nLoci){
		localVector<-outputVector[which(outputIndex==i)]
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
	lnL.byLocus<-unname(meanIntercepts + (totalPopVector * meanSlopes))
	names(lnL.byLocus)<-paste("lnL.L",c(1:length(lnL.byLocus)),sep="")
	lnL<-sum(lnL.byLocus)
	names(lnL)<-"lnL"
	#Prep to return
	if(saveNoExtrap==TRUE){
		lnL<-c(lnL,lnL.noExtrap=lnL.noExtrap,lnL.byLocus,meanSlopes,meanIntercepts,coeffs95Lower,coeffs95Upper)
	}else{
		lnL<-c(lnL,lnL.byLocus,meanSlopes,meanIntercepts,coeffs95Lower,coeffs95Upper)
	}
	lnL.mat<-data.frame(matrix(nrow=1,ncol=length(lnL),lnL))
	colnames(lnL.mat)<-names(lnL)

	return(lnL.mat)
}

#If only one set of subsamples is used, use this to calculate likelihoods (no extrapolation),
#although effective subsample sizes can be considered by using SumDivScaledNreps
ConvertOutputVectorToLikelihood.1sub<-function(outputVector,popAssignments,nTrees=1,subsamplesPerGene=1, 		
		totalPopVector,summaryFn="mean", nEq=100) {
	outputVector<-as.numeric(outputVector)
	outputVector<-AdjustUsingBeta(numMatches=outputVector, nTrees=nTrees, nTips=sum(popAssignments[[1]]), nEq=nEq)
	outputVector<-log(outputVector)
	lnL<-0
	localVector<-rep(NA, subsamplesPerGene)
	baseIndex<-1
	locusIndex<-1
	for (i in sequence(length(outputVector))) {
		localVector[baseIndex]<-outputVector[i]
#		print(localVector)
		baseIndex <- baseIndex+1
		if(i%%subsamplesPerGene == 0) {
			if(summaryFn=="SumDivScaledNreps"){
				lnL<-lnL+get(summaryFn)(localVector=localVector,popAssignments,totalPopVector=totalPopVector[locusIndex],
					subsamplesPerGene=subsamplesPerGene)
			}else{
				lnL<-lnL+get(summaryFn)(localVector)
			}
			localVector<-rep(NA, subsamplesPerGene)
			baseIndex<-1
			locusIndex<-locusIndex + 1
		}
	}
	return(lnL)
}





#GetAndBindLabel and GetClades together get a list of all clades in the tree. For the taxa descended from each clade, 
#sorts alphabetically and then makes them a string. These are slow, and have since been replaced by GetCladesQuickly
GetClades <- function(phy) {
	return(simplify2array(sapply(subtrees(phy), GetAndBindLabel)))
}

#Used with GetClades
GetAndBindLabel <- function(phy) { 
	#note the sorting here
	return( paste( sort( phy$tip.label ), collapse="_" ) )
}

#This extracts clades from raw newick trees (implemented for efficiency). It also allows for
#relabeling of individuals with population letters. Thus, this replaces GetCladesQuickly (7 times
#slower) and ConvertAlleleNumbersToPopulationLetters
GetCladeNamesFromNewickString <- function(newick.string,assignFrame,convertIndivNumsToPopLetters=FALSE,
		getOutDegreeOfPolytomies=FALSE) {
	descendantCounts<-c() #for saving the degree of ploytomies
	nobrlen.tree <- gsub(pattern=";", replacement="", x=gsub(pattern='\\:\\d+\\.*\\d*', replacement="", 
		x= newick.string, perl=TRUE)) #remove branch lengths
	clade.names<-rep(NA, sum(grepl("\\(", strsplit(nobrlen.tree,"")[[1]]))) #num of clades=num of parentheses/2
	
	#Start the culling of clades
	for(clade.count in sequence(length(clade.names))) {
		m<-regexpr(pattern="\\([^(^)]+\\)", text=nobrlen.tree, perl=TRUE) #find location of first set of closed parentheses
		clade<-gsub("\\)", "", gsub("\\(","",regmatches(x=nobrlen.tree, m))) #everything in between is a clade
		taxa<-strsplit(clade, ",")[[1]]

		#Get the degree of ploytomies for each clade with more than two taxa
		if(getOutDegreeOfPolytomies){
			if(length(taxa) > 2){ #if there are more than 2 members a new clade...
				descendantCounts<-append(descendantCounts,length(taxa)) #record the degree
			}
		}
		
		#Convert numbers to letters
		if(convertIndivNumsToPopLetters){
			which.untransformed.alleles<-which(grepl("\\d", taxa, perl=TRUE)) #number of untransformed
			for(allele.index in sequence(length(which.untransformed.alleles))) {
				taxon.index<-which.untransformed.alleles[allele.index]
				taxa[taxon.index]<-as.character(assignFrame$popLabel[match(taxa[taxon.index],assignFrame$indivTotal)])
			}
		}
		
		clade.name<-paste(sort(taxa),collapse="-")
		regmatches(x=nobrlen.tree, m)<-clade.name
		clade.names[clade.count]<-clade.name	
	}
	if(getOutDegreeOfPolytomies){
		return(list(clade.names,descendantCounts))
	}else{
		return(clade.names)
	}
}	

GetOutDegreeOfPolytomies <- function(phy) {
	descendantCounts <- table(phy$edge[,1])
	descendantCounts <- unname(descendantCounts[which(descendantCounts>2)])
	return(descendantCounts)
}

#Match clades of two trees
#assumes:
#1. You have already run ConvertAlleleNumbersToPopulationLetters so these trees have letters
#2. You have already made them have the same size (do SubsampleMSTree if needed)
GetScoreOfSingleTree <- function(cladesMS, phyMS, cladesGene, phyGene, polytomyDegreesForCorrection=NULL) {
	numberCladesInMSOnly <- sum(!cladesMS%in%cladesGene)
	numberCladesInGeneOnly <- sum(!cladesGene%in%cladesMS)
	matchCount <- 0
	if(numberCladesInMSOnly == 0 && numberCladesInGeneOnly==0) {
		matchCount <- 1
	}
	if(!is.null(polytomyDegreesForCorrection)){
		if(numberCladesInMSOnly > 0 && numberCladesInGeneOnly == 0){
			descendantCounts<-polytomyDegreesForCorrection
			correction <- 1
			for (i in sequence(length(descendantCounts))) {
				correction <- correction * howmanytrees(descendantCounts[i], rooted=TRUE, binary=TRUE, labeled=TRUE)
			}
			matchCount <- 1 / correction #idea here is that the gene tree could have been resolved at each polytomy multiple ways
			#only one of these ways would match the given phyMS tree. So we figure out the number of ways to resolve polytomy 1, 
			#multiply that by the number of ways to resolve polytomy 2, etc. A polytomy with three descendant edges has 3 ways to 
			#resolve it, one with 4 descendant edges has 3 * 5 = 15, etc.
		}
	}
	return(matchCount)
}

#Drops tips from the simulated trees (to be used in accordance with popAssignments)
#Relies on using phylo objects, which makes it slow
SubsampleMSTree<-function(phy,popVectorOrig,popVectorFinal) {
	taxaToDrop<-c()
	minSample<-1
	maxSample<-0
	uniquePops<-unique(phy$tip.label)
	for(population in 1:length(uniquePops)){		
		maxSample<-minSample + popVectorOrig[population] - 1
		sampleSizeDiff<-popVectorOrig[population] - popVectorFinal[population]
		if(sampleSizeDiff > 0){
			taxaToDrop<-append(taxaToDrop, 
				which(phy$tip.label==uniquePops[population])[c(1:sampleSizeDiff)])
		}
		minSample<-maxSample + 1
	}
	phy<-drop.tip(phy,taxaToDrop)

	return(phy)
}

#Drop tips from newick string. This is not presently working. If we just keep all the 
#parentheses as is when deleting taxa and just drop labels and commas, this will be easy to get working.
SubsampleMSTree.raw<-function(phy,popVectorOrig,popVectorFinal) {
	#Get taxon labels from tree
	taxonLabels<-strsplit(phy, "[^\\w]", perl=TRUE)[[1]]
	taxonLabels<-taxonLabels[which(nchar(taxonLabels)>0)]

	#Get tree positions for labels to drop
	uniquePops<-unique(taxonLabels)
	
	#Split tree into characters and match to old taxon label positions
	split.tree<-split.tree<-strsplit(phy,"")[[1]]
	labelPositions<-c()
	for(i in 1:length(uniquePops)){
		labelPositions<-append(labelPositions,which(split.tree==uniquePops[i]))
	}
	
	#Select positions in the split tree of taxa to drop
	taxaToDrop<-c()
	startPos<-1
	endPos<-popVectorOrig[population]
	for(population in 1:length(uniquePops)){		
		sampleSizeDiff<-popVectorOrig[population] - popVectorFinal[population]
		if(sampleSizeDiff > 0){
			currentLabels<-labelPositions[startPos:endPos]
			taxaToDrop<-append(taxaToDrop,currentLabels[c(sample(1:popVectorOrig[population],sampleSizeDiff,replace=FALSE))])
			startPos<-startPos + popVectorOrig[population]
			endPos<-endPos + popVectorOrig[population]
		}
	}

	#Drop one tip at a time (must be fully-resolved)	
	#Toss the label and adjacent commas and opening brackets
	for(h in 1:length(taxaToDrop)){
		split.tree[taxaToDrop[h]]<-"drop"
		if(split.tree[taxaToDrop[h] - 1] == "("){
			split.tree[taxaToDrop[h] - 1]<-"drop"
			bracket="left"
		}
		if(split.tree[taxaToDrop[h] - 1] == ","){
			split.tree[taxaToDrop[h] - 1]<-"drop"
		}
		if(split.tree[taxaToDrop[h] + 1] == ")"){
			split.tree[taxaToDrop[h] + 1]<-"drop"
			bracket="right"
		}
		if(split.tree[taxaToDrop[h] + 1] == ","){
			split.tree[taxaToDrop[h] + 1]<-"drop"
		}
		
		#Find and toss the closing bracket
		if(bracket=="right"){
			leftCount<- 0 #when leftCount exceeds rightCount, drop parenthesis
			rightCount<- 0
			currentRemainder<-split.tree[1:(taxaToDrop[h] - 2)]
			rightPos<-which(currentRemainder==")")
			leftPos<-which(currentRemainder=="(")
			continue<-TRUE
			for(j in length(currentRemainder):1){
				if(continue==TRUE){
					if(j%in%leftPos){
						leftCount<-leftCount + 1
					}
					if(j%in%rightPos){
						rightCount<-rightCount + 1
					}
					if(leftCount > rightCount){
						split.tree[j]<-"drop"
						continue<-FALSE
					}
				}
			}
		}else{
			leftCount<-0 #when rightCount exceeds leftCount, drop parenthesis
			rightCount<-0
			currentRemainder<-split.tree[(taxaToDrop[h] + 2):length(split.tree)]
			positionRemainder<-length(split.tree[1:(taxaToDrop[h] + 2)])
			rightPos<-which(currentRemainder==")")
			leftPos<-which(currentRemainder=="(")
			continue<-TRUE
			for(j in 1:length(currentRemainder)){
				if(continue==TRUE){
					if(j%in%leftPos){
						leftCount<-leftCount + 1
					}
					if(j%in%rightPos){
						rightCount<-rightCount + 1
					}
					if(rightCount > leftCount){
						split.tree[j + (positionRemainder - 1)]<-"drop"
						continue<-FALSE
					}
				}
			}
		}
	}

	#Drop tips
	split.tree<-split.tree[which(split.tree != "drop")]

	#Get rid of extra commas
	if(split.tree[1]==","){
		split.tree[1]<-"drop"
	}
	if(split.tree[length(split.tree)]==","){
		split.tree[length(split.tree)]<-"drop"
	}
	for(l in 1:2){ #twice to be sure
		for(k in 1:(length(split.tree) - 1)){
			if(split.tree[k]=="(" && split.tree[k+1]==","){
				split.tree[k+1]<-"drop"
			}
			if(split.tree[k]=="," && split.tree[k+1]==")"){
				split.tree[k]<-"drop"
			}
			if(split.tree[k]=="," && split.tree[k+1]==","){
				split.tree[k]<-"drop"
			}
		}
		split.tree<-split.tree[which(split.tree != "drop")]
	}
	split.tree<-split.tree[which(split.tree != "drop")]
	split.tree.cat<-paste(split.tree,collapse="")


	save<-split.tree.cat #for troubleshooting
	split.tree.cat<-save #for troubleshooting

	#Get rid of brackets around singletons
	while(length(grep("\\(\\(\\w\\)",split.tree.cat)) == 1){
		m<-regexpr(pattern="\\(\\(\\w\\)",text=split.tree.cat, perl=TRUE)
		new.clade<-gsub("\\)","",regmatches(x=split.tree.cat, m))
		if((m[1] + 3) < length(split.tree)){
			split.tree<-c(split.tree[1:(m[1] - 1)],strsplit(new.clade,"")[[1]],
				split.tree[(m[1] + 4):length(split.tree)])
		}else{
			split.tree<-c(split.tree[1:(m[1] - 1)],strsplit(new.clade,"")[[1]])
		}
		split.tree.cat<-paste(split.tree,collapse="")
	}

	while(length(grep("\\(\\w\\)\\)",split.tree.cat)) == 1){
		m<-regexpr(pattern="\\(\\w\\)\\)",text=split.tree.cat,perl=TRUE)
		new.clade<-gsub("\\(","",regmatches(x=split.tree.cat,m))
		split.tree<-strsplit(split.tree.cat,"")[[1]]
		if((m[1]) > 1){
			split.tree<-c(split.tree[1:(m[1] - 1)],strsplit(new.clade,"")[[1]],
				split.tree[(m[1] + 4):length(split.tree)])
		}else{
			split.tree<-c(strsplit(new.clade,"")[[1]],split.tree[2:length(split.tree)])
		}
		split.tree.cat<-paste(split.tree,collapse="")
	}

	while(length(grep("\\(\\w\\)",split.tree.cat)) == 1){
		m<-regexpr(pattern="\\(\\w\\)",text=split.tree.cat,perl=TRUE)[1]
		split.tree.temp<-strsplit(split.tree.cat,"")[[1]]
		opening.l<-length(split.tree.temp[1:(m-1)][which(split.tree.temp=="(")])
		closing.l<-length(split.tree.temp[(m+3):length(split.tree.temp)][which(split.tree.temp==")")])
		if(opening.l > closing.l || m==length(split.tree)){
			split.tree.temp<-split.tree.temp[-m]
		}
		if(closing.l > opening.l || m==1){
			split.tree.temp<-split.tree.temp[-(m+2)]
		}
		split.tree.cat<-paste(split.tree.temp,collapse="")
	}	
	
	#Get rid of double brackets	
	#First get positions of double brackets
	openBracket<-c()
	closeBracket<-c()
	split.tree<-strsplit(split.tree.cat,"")[[1]]
	for(m in 1:(length(split.tree) - 1)){
		if(split.tree[m]=="(" && split.tree[m+1]=="("){
			openBracket<-append(openBracket,m)
		}
		if(split.tree[m]==")" && split.tree[m+1]==")"){
			closeBracket<-append(closeBracket,m + 1)
		}
	}
	
	#For different pairs, if there are as many or more brackets as taxa, toss them
	end<-FALSE
	while((length(grep("\\(",split.tree)) >= length(grep("\\w",split.tree)) || 
			length(grep("\\)",split.tree)) >= length(grep("\\w",split.tree))) &&
			end==FALSE){	
		if(openBracket[length(openBracket)] > closeBracket[1]){
			substring<-split.tree[openBracket[1]:closeBracket[length(closeBracket)]]
			if(length(grep("\\(",split.tree)) >= length(grep("\\w",split.tree))){
				if(length(grep("\\(",substring)) >= length(grep("\\w",substring))){
					substring<-substring[-1]
					removed<-TRUE
				}
			}
			
			if(length(grep("\\)",split.tree)) >= length(grep("\\w",split.tree))){
				if(length(grep("\\)",substring)) >= length(grep("\\w",substring))){
					substring<-substring[-length(substring)]
					removed<-TRUE
				}
			}
			if(removed==TRUE){		
				split.tree<-split.tree[-c(openBracket[1]:closeBracket[length(closeBracket)])]
				if(openBracket[1] > 1){
					split.tree<-c(split.tree[1:openBracket[1] - 1],substring,
						split.tree[openBracket[1]:length(split.tree)])
				}else{
					split.tree<-c(substring,split.tree)
				}
				if(length(openBracket) > 1 && length(closeBracket) > 1){
					openBracket<-openBracket[-1]
					closeBracket<-closeBracket[-(length(openBracket))]	
				}else{
					end<-TRUE
				}
			}else{
				end<-TRUE
			}
			
		}else{
			removed<-FALSE
			substring<-split.tree[openBracket[length(openBracket)]:closeBracket[1]]
			if(length(grep("\\(",split.tree)) >= length(grep("\\w",split.tree))){
				if(length(grep("\\(",substring)) >= length(grep("\\w",substring))){
					substring<-substring[-1]
					removed<-TRUE
				}
			}
			if(length(grep("\\)",split.tree)) >= length(grep("\\w",split.tree))){
				if(length(grep("\\)",substring)) >= length(grep("\\w",substring))){
					substring<-substring[-length(substring)]
					removed<-TRUE
					reduce<-TRUE
				}
			}
			
			if(removed==TRUE){	
				split.tree<-split.tree[-c(openBracket[length(openBracket)]:closeBracket[1])]
				if(openBracket[length(openBracket)] > 1){
					split.tree<-c(split.tree[1:openBracket[length(openBracket)] - 1],substring,
					split.tree[openBracket[length(openBracket)]:length(split.tree)])
				}else{
					split.tree<-c(substring,split.tree)
				}
				if(length(openBracket) > 1 && length(closeBracket) > 1){
					openBracket<-openBracket[-(length(openBracket))]
					closeBracket<-closeBracket[-1]
					if(reduce==TRUE){
						closeBracket<-closeBracket - 1
						reduce==FALSE
					}
				}else{
					end<-TRUE
				}
			}else{
				end<-TRUE
			}
		}
	}	
	split.tree.cat<-paste(split.tree,collapse="")
	
	return(split.tree.cat)
}