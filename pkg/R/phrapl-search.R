#TO DO: If the optimal value is outside the bounds of the grid, offer warning or option to restart search centered at new grid
#I have tweaked the code for SearchContinuousSpaceNLopter such that 1) different grid parameters were used and such that 2)the grid file is saved
#with a different name for each run
SearchContinuousModelSpaceNLoptr<-function(p, migrationArrayMap, migrationArray, popAssignments, badAIC=100000000000000, 
	maxParameterValue=100, nTrees=2e5, nTreesGrid=nTrees ,msPath="ms",comparePath=system.file("extdata", "comparecladespipe.pl", package="phrapl"),
	subsampleWeights.df=NULL,
	unresolvedTest=TRUE,print.ms.string=FALSE, ncores=1,print.results=FALSE,print.matches=FALSE,debug=FALSE,method="nlminb",
	itnmax=NULL, return.all=FALSE, maxtime=0, maxeval=0, parameterBounds=list(minCollapseTime=0.1, minCollapseRatio=0, 
	minN0Ratio=0.1, minMigrationRate=0.05, minMigrationRatio=0.1), numReps=0, startGrid=startGrid, 
	collapseStarts=c(0.30,0.58,1.11,2.12,4.07,7.81,15.00), n0Starts=c(0.1,0.5,1,2,4), 
	migrationStarts=c(0.10,0.22,0.46,1.00,2.15,4.64,10.00), gridSave=NULL,gridSaveFile=NULL,subsamplesPerGene=1,
	totalPopVector,summaryFn="mean",saveNoExtrap=FALSE, doSNPs=FALSE, nEq=100, ...) {
  modelID<-ReturnModel(p,migrationArrayMap)
  best.result <- c()
  best.result.objective <- badAIC
  if(print.results) {
#    resultVector<-c(modelID,p)
#    names(resultVector)<-c("migrationArryIndividualID","collapseMatrix.number", "n0multiplierMap.number","migrationArray.number")
#    print(resultVector)
  }
  if(is.na(modelID)) {
    return(badAIC)
  }
  else {
    paramNames<-MsIndividualParameters(migrationArray[[modelID]])
    if(is.null(startGrid)) { #need to make our own grid
    	startingVectorList<-list()
    	nameCount<-1
    	numCollapse <- sum(grepl("collapse",paramNames))
    	if (numCollapse > 0) {
    		for (i in sequence(numCollapse)) {
    			startingVectorList<-append(startingVectorList, list(collapseStarts))
    			names(startingVectorList)[nameCount]<-grep("collapse",paramNames,value=TRUE)[i]
    			nameCount<-nameCount + 1
    		} 
    	}
        numn0 <- sum(grepl("n0multiplier",paramNames))-1 #-1 since one is fixed to be 1
    	if (numn0 > 0) {
    		for (i in sequence(numn0)) {
    			startingVectorList<-append(startingVectorList, list(n0Starts))
    			names(startingVectorList)[nameCount]<-grep("n0multiplier",paramNames,value=TRUE)[i]
    			nameCount<-nameCount + 1
    		} 
    	}	
    	numMigration <- sum(grepl("migration",paramNames))
    	if (numMigration > 0) {
    		for (i in sequence(numMigration)) {
    			startingVectorList<-append(startingVectorList, list(migrationStarts))
    			names(startingVectorList)[nameCount]<-grep("migration",paramNames,value=TRUE)[i]
    			nameCount<-nameCount + 1
    		} 
    	}
    	startGrid <- CreateStartGrid(startingVectorList)
    	startGrid<-startGrid[[1]] #default grid shouldn't be list (as there is always one grid)
    }

    if(is.null(nTreesGrid)) {
    	nTreesGrid<-10*nTrees #thinking here that want better estimate on the grid than in the heat of the search
    }
  	
  	#Get and store AIC for each set of grid values
  	initial.AIC<-c()
  	for(t in 1:nrow(startGrid)){
    	currentAIC<-ReturnAIC(par=as.numeric(startGrid[t,]),migrationIndividual=migrationArray[[modelID]],
    	badAIC=badAIC,maxParameterValue=maxParameterValue,nTrees=nTreesGrid,msPath=msPath,comparePath=comparePath,
    	unresolvedTest=unresolvedTest,print.ms.string=print.ms.string,print.results=print.results,print.matches=print.matches,
    	ncores=ncores,debug=debug,numReps=numReps,parameterBounds=parameterBounds,subsamplesPerGene=subsamplesPerGene,
    	totalPopVector=totalPopVector,popAssignments=popAssignments,subsampleWeights.df=subsampleWeights.df,
    	summaryFn=summaryFn,saveNoExtrap=saveNoExtrap,doSNPs=doSNPs,nEq=nEq)
    	initial.AIC<-append(initial.AIC,currentAIC)
    }
 
    if(debug) {
    	print(cbind(AIC=initial.AIC,exp(startGrid)))
    }

	#Save the grid as object
    positionOfFirstN0 <- min(grep("n0multiplier", MsIndividualParameters(migrationArray[[modelID]])))
    if(length(popAssignments) > 1){
   	 	thisGrid<-cbind(AIC=initial.AIC[,1],exp(startGrid),initial.AIC[,2:length(initial.AIC)])
   	}else{
   	 	thisGrid<-cbind(AIC=initial.AIC,exp(startGrid))
   	}
    thisGrid<-thisGrid[order(thisGrid$AIC),]
    #Save grid to file
    if(!is.null(gridSave)) {
    	save(thisGrid, file=gridSave)
    }
    for(rep in sequence(min(numReps, dim(startGrid)[1]))) {
 	   #startingVals<-log(c(rlnorm(sum(grepl("collapse",paramNames)),1,1), rlnorm(-1+sum(grepl("n0multiplier",paramNames)),1,1), rbeta(sum(grepl("migration",paramNames)),shape1=1,shape2=3) )) #going to optimize in log space. Remove the first n0 parameter from optimization vector
 	   startingVals <- unlist(startGrid[order(initial.AIC)[rep],]) #order(initial.AIC) gives indices of best values, min to max, so if we want the second best it's the second one here. NA are stuck last, if present
 	   if(debug) {
 	     print(startingVals) 
    	}
		#Currently, nloptr optimization will not work because ReturnAIC is returning a list that includes the AIC plus summary information about the slopes used
		#to estimate AIC (the eval_f can only use a function that returns AIC alone). If we want to save slope information while using nloptr, we could fix 
		#ReturnAIC to return just AIC, but also to store the best slopes list using a global variable (something like "if(currentAIC[[1]]<previousAIC[[1]]){
		#previousAIC<<-currentAIC}, which could be saved
 	   	searchResults<-nloptr(x0=startingVals, eval_f=ReturnAIC, opts=list("maxeval"=itnmax, "algorithm"="NLOPT_LN_SBPLX", "print_level"=1,
 	   		maxtime=maxtime, maxeval=maxeval), migrationIndividual=migrationArray[[modelID]], popVector=popVector, 
 	   		badAIC=badAIC, maxParameterValue=maxParameterValue,
 	   		nTrees=nTrees,msPath=msPath,comparePath=comparePath,unresolvedTest=unresolvedTest,ncores=ncores,
 	   		print.ms.string=print.ms.string, print.results=print.results,print.matches=print.matches,
 	   		debug=debug,numReps=numReps,parameterBounds=parameterBounds,subsamplesPerGene=subsamplesPerGene,summaryFn=summaryFn,
 	   		totalPopVector=totalPopVector,subsampleWeights.df=subsampleWeights.df,popAssignments=popAssignments,
 	   		saveNoExtrap=saveNoExtrap, doSNPs=doSNPs, nEq=nEq)
 	  

  		#stitch the first N0multiplier (=1) into the final parameter vector
		positionOfFirstN0 <- min(grep("n0multiplier", MsIndividualParameters(migrationArray[[modelID]])))
		solutionVectorFirstPart<-searchResults$solution[sequence(positionOfFirstN0-1)]
		solutionVectorSecondPart<-searchResults$solution[(1+length(solutionVectorFirstPart)):length(searchResults$solution)]
		if((1+length(solutionVectorFirstPart)) > length(searchResults$solution)) {
  			solutionVectorSecondPart<-c()
  		}
		searchResults$solution<-c(solutionVectorFirstPart,log(1),solutionVectorSecondPart)
		 if(debug) {
# 	   	print(searchResults)
 	   	}
 	   	if(searchResults$objective <= best.result.objective) {
 	   		best.result<-searchResults
 	   		best.result.objective<-searchResults$objective	
 	   }
    }
	best.resultGrid<-list(best.result,thisGrid)
    ifelse(return.all, return(best.resultGrid), return(best.result.objective))     
  }
}

#This function was formally known as "ExhaustiveSearchNLoptr", and still has the capabilities of running a heuristic
#nloptr search. However, since we are currently focusing on a grid search, we have changed the name of the function
GridSearch<-function(modelRange=c(1:length(migrationArray)), migrationArrayMap,migrationArray,popAssignments, 
		badAIC=100000000000000,maxParameterValue=100,nTrees=1e5 ,msPath="ms",comparePath=system.file("extdata", "comparecladespipe.pl", package="phrapl"),
		observedTrees=NULL,subsampleWeights.df=NULL, doWeights=TRUE, unresolvedTest=TRUE,
		print.ms.string=FALSE,print.results=FALSE,print.matches=FALSE,debug=FALSE,method="nlminb",itnmax=NULL,
		ncores=1,results.file=NULL,maxtime=0, maxeval=0,return.all=TRUE, numReps=0,startGrid=NULL,
		collapseStarts=c(0.30,0.58,1.11,2.12,4.07,7.81,15.00), n0Starts=c(0.1,0.5,1,2,4), 
		migrationStarts=c(0.10,0.22,0.46,1.00,2.15,4.64,10.00),subsamplesPerGene=1,
		totalPopVector=NULL,summaryFn="mean",saveNoExtrap=FALSE, doSNPs=FALSE, nEq=100, ...){

	if(length(modelRange) != length(migrationArray)) { #need to look at only particular rows
		migrationArray<-migrationArray[modelRange]
		migrationArrayMap<-GenerateMigrationArrayMapTrunc(migrationArrayMap,modelRange)
	}
	if(is.null(subsampleWeights.df) && doWeights) {
   	subsampleWeights.df <- GetPermutationWeightsAcrossSubsamples(popAssignments, observedTrees)
	}
  	#Prepare temporary tree and assign files
	for(k in 1:length(popAssignments)){
		write.table(CreateAssignment.df(popAssignments[[k]]),file=paste(tempdir(),"/assign",k,".txt",sep=""),
			quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=FALSE)
		if(file.exists(paste(tempdir(),"/observed",k,".tre",sep=""))){
			unlink(paste(tempdir(),"/observed",k,".tre",sep=""))
		}
	}
	for(f in 1:length(observedTrees)){ #for each popAssignment
		if (class(observedTrees[[f]]) != "multiPhylo"){
        	observedTrees[[f]] <- c(observedTrees[[f]])
    	}
		for(m in 1:length(observedTrees[[f]])){ #for each tree
			write.tree(observedTrees[[f]][[m]],file=paste(tempdir(),"/observed",f,".tre",sep=""),append=TRUE)
		}		
	}
 
  AIC.values<-rep(NA,length(migrationArray))
  gridList<-list() #for storing model grids
  results.list<-list()
  for (i in sequence(length(migrationArray))) {
  	p<-c(migrationArrayMap$collapseMatrix.number[i], migrationArrayMap$n0multiplierMap.number[i], migrationArrayMap$migrationArray.number[i])
  	if(return.all) {
 		if(!is.null(startGrid)){
 			currentStartGrid<-startGrid[[i]] #imported startGrid must be a list
 		}else{
 			currentStartGrid=NULL
 		}

  		result.indiv<-NULL
  		try(result.indiv<-SearchContinuousModelSpaceNLoptr(p,migrationArrayMap,migrationArray,popAssignments,badAIC=badAIC,
  		maxParameterValue=maxParameterValue,nTrees=nTrees,msPath=msPath,comparePath=comparePath,unresolvedTest=unresolvedTest,
  		print.ms.string=print.ms.string,print.results=print.results,print.matches=print.matches,debug=debug,
  		method=method,itnmax=itnmax, maxtime=maxtime,ncores=ncores,maxeval=maxeval, return.all=return.all,numReps=numReps,
  		startGrid=currentStartGrid,gridSave=NULL,collapseStarts=collapseStarts,n0Starts=n0Starts,migrationStarts=migrationStarts,
  		subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,subsampleWeights.df=subsampleWeights.df,
  		summaryFn=summaryFn,saveNoExtrap=saveNoExtrap, doSNPs=doSNPs, nEq=nEq,...))
  		
		gridList[[length(gridList)+1]]<-result.indiv[[2]] #make list of model grids
#  		print(result.indiv[[1]])
  		if(!is.null(result.indiv[[1]])) {
  			AIC.values[i]<-result.indiv[[1]]$objective	
  		}
  		results.list<-append(results.list, list(result.indiv[[1]]))
  	} else {
  		try(AIC.values[i]<-SearchContinuousModelSpaceNLoptr(p,migrationArrayMap,migrationArray,popAssignments,badAIC=badAIC,
  		maxParameterValue=maxParameterValue, nTrees=nTrees, msPath=msPath,comparePath=comparePath,ncores=ncores,
  		unresolvedTest=unresolvedTest,print.ms.string=print.ms.string,print.results=print.results,print.matches=print.matches,
  		debug=debug,method=method,itnmax=itnmax, maxtime=maxtime, maxeval=maxeval, return.all=return.all,numReps=numReps,
  		startGrid=currentStartGrid,collapseStarts=collapseStarts,n0Starts=n0Starts,migrationStarts=migrationStarts,
  		gridSave=NULL,subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,subsampleWeights.df=subsampleWeights.df,
  		summaryFn=summaryFn,saveNoExtrap=saveNoExtrap, doSNPs=doSNPs, nEq=nEq, ...))
  	}

#  	 	print(c(i, length(migrationArray), i/length(migrationArray), AIC.values[i]))

  		if(!is.null(results.file)) {
			save(list=ls(), file=results.file)
		}
	}

  	#Toss temporary tree and assign files
	for(k in 1:length(popAssignments)){
		unlink(c(paste(tempdir(),"/assign",k,".txt",sep=""),paste(tempdir(),"/observed",k,".tre",sep=""),paste(tempdir(), "/mstrees.txt", sep="")))
	}

	#Save table of best models
	if(numReps==0){
		overall.results<-ExtractGridAICs(result=gridList,migrationArray=migrationArray)
	}else{
		overall.results<-ExtractAICs(result=results.list,migrationArray=migrationArray)
	}
	overall.results$models <- modelRange #so if we've done subsampling, use right model numbers
	#Save parameter estimates and parameter indexes to tables
	if(numReps==0){
		parameters<-ExtractGridParameters(migrationArray=migrationArray,result=gridList,
			popVector=popAssignments[[1]])
	}else{
		parameters<-ExtractParameters(migrationArray=migrationArray,result=results.list,
			popVector=popAssignments[[1]])
	}
	
	#Add all results to list
	if(numReps==0){
		results.final<-list("AIC.Grid"=gridList,"overall.results"=overall.results,
			"parameters"=parameters[[1]],"parameterIndexes"=parameters[[2]])
	}else{
		results.final<-list("search.results"=results.list,"AIC.Grid"=gridList,
			"overall.results"=overall.results,"parameters"=parameters[[1]],"parameterIndexes"=parameters[[2]])
	}
	
	ifelse(return.all, return(results.final), return(AIC.values))     
}