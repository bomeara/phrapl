#TO DO: If the optimal value is outside the bounds of the grid, offer warning or option to restart search centered at new grid
SearchContinuousModelSpaceNLoptr<-function(p, migrationArrayMap=NULL, migrationArray, popAssignments, badAIC=100000000000000, 
	maxParameterValue=20, nTrees=2e5, nTreesGrid=nTrees ,msPath="ms",comparePath=system.file("extdata", "comparecladespipe.pl", package="phrapl"),
	subsampleWeights.df=NULL,
	unresolvedTest=TRUE,print.ms.string=FALSE, ncores=1,print.results=TRUE,print.matches=FALSE,debug=FALSE,method="nlminb",
	itnmax=NULL, return.all=FALSE, maxtime=0, maxeval=0, parameterBounds=list(minCollapseTime=0.1,
	minCollapseRatio=0,minN0Ratio=0.1,minGrowth=0.1,minGrowthRatio=0.1,minMigrationRate=0.05,minMigrationRatio=0.1), 
	numReps=0, startGrid=startGrid, collapseStarts=c(0.30,0.58,1.11,2.12,4.07,7.81,15.00), n0Starts=c(0.1,0.5,1,2,4), 
	growthStarts=c(0.30,0.58,1.11,2.12,4.07,7.81,15.00),migrationStarts=c(0.10,0.22,0.46,1.00,2.15,4.64), gridSave=NULL,
	gridSaveFile=NULL,subsamplesPerGene=1,totalPopVector,summaryFn="mean",saveNoExtrap=FALSE,doSNPs=FALSE,nEq=100,
	setCollapseZero=NULL,rm.n0=TRUE,popScaling=NULL,checkpointFile=NULL,addedEventTime=NULL,addedEventTimeAsScalar=TRUE,
	optimization="grid",...) {
	if(!is.null(migrationArrayMap)){
		modelID<-ReturnModel(p,migrationArrayMap)
  	}else{
  		modelID<-p
	}
	best.result <- c()
	best.result.objective <- badAIC
	if(print.results) {
#    resultVector<-c(modelID,p)
#    names(resultVector)<-c("migrationArryIndividualID","collapseMatrix.number", "n0multiplierMap.number","migrationArray.number")
#    print(resultVector)
	}
	if(is.na(modelID)) {
    	return(badAIC)
	}else{
    	paramNames<-MsIndividualParameters(migrationArray[[modelID]])
    	if(is.null(startGrid)) { #need to make our own grid
    		startingVectorList<-list()
    		nameCount<-1
    	
			#Calculate number of collapse parameters to estimate (subtracting off those that are set to zero)
			numCollapse <- sum(grepl("collapse",paramNames)) - length(setCollapseZero)
			if(numCollapse < 0){ #If for some reason there are more fixed parameters than possible, reduce them
				setCollapseZero<-sequence(sum(grepl("collapse",paramNames)))
				numCollapse <- sum(grepl("collapse",paramNames)) - length(setCollapseZero)
			}
			if(!is.null(setCollapseZero)){ #If setCollapseZero is nonzero
				collapses2estimate<-grep("collapse",paramNames,value=TRUE)[which(!grep("collapse",paramNames) %in% setCollapseZero)]
			}else{
				collapses2estimate<-grep("collapse",paramNames,value=TRUE)
			}
			if(numCollapse > 0) {
				for (i in sequence(numCollapse)) {
					startingVectorList<-append(startingVectorList, list(collapseStarts))
					names(startingVectorList)[nameCount]<-collapses2estimate[i]
					nameCount<-nameCount + 1
				} 
			}
    	
			#Calculate number of n0 parameters (if just one, subtract it off)
			if(sum(grepl("n0multiplier",paramNames)) == 1){
			   numn0 <- sum(grepl("n0multiplier",paramNames)) - 1 # -1 since one is fixed to be 1
			}else{
			   numn0 <- sum(grepl("n0multiplier",paramNames))
			}
			if (numn0 > 0) {
				for (i in sequence(numn0)) {
					startingVectorList<-append(startingVectorList, list(n0Starts))
					names(startingVectorList)[nameCount]<-grep("n0multiplier",paramNames,value=TRUE)[i]
					nameCount<-nameCount + 1
				} 
			}	
    	
			#Calculate number of growth parameters
			numGrowth <- sum(grepl("growth",paramNames))
			if (numGrowth > 0) {
				for (i in sequence(numGrowth)) {
					startingVectorList<-append(startingVectorList, list(growthStarts))
					names(startingVectorList)[nameCount]<-grep("growth",paramNames,value=TRUE)[i]
					nameCount<-nameCount + 1
				} 
			}
    	
			#Calculate number of migration parameters
			numMigration <- sum(grepl("migration",paramNames))
			if (numMigration > 0) {
				for (i in sequence(numMigration)) {
					startingVectorList<-append(startingVectorList, list(migrationStarts))
					names(startingVectorList)[nameCount]<-grep("migration",paramNames,value=TRUE)[i]
					nameCount<-nameCount + 1
				} 
			}
    	
			#If some collapses are set to zero, there are collapses in the model, and also added events,
			#wait to filter the grid for added events until after zero-fixed collapses are stitched in below
			if(!is.null(setCollapseZero) && sum(grepl("collapse",paramNames)) != 0 && !is.null(addedEventTime)){
				startGrid<-CreateStartGrid(fineGrid=startingVectorList,migrationIndividual=migrationArray[[modelID]])
				startGrid<-startGrid[[1]] #default grid shouldn't be list (as there is always one grid)
		
			#Otherwise, filter for any added events now
			}else{
				startGrid<-CreateStartGrid(fineGrid=startingVectorList,migrationIndividual=migrationArray[[modelID]],
					addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar)
				startGrid<-startGrid[[1]] #default grid shouldn't be list (as there is always one grid)
			}
    	
			#If some collapses are set to zero (and there are collapses in the model), stitch these values into the startGrid
			if(!is.null(setCollapseZero) && sum(grepl("collapse",paramNames)) != 0){
			
				#If there are only fixed parameters...
				if(length(startGrid) == 0){
					startGrid<-data.frame(matrix(rep(log(0),length(setCollapseZero)),nrow=1,ncol=length(setCollapseZero)))
					for(g in 1:length(setCollapseZero)){
						names(startGrid)[g]<-paste("collapse_",setCollapseZero[g],sep="")
					}
				}else{
			
					#Otherwise, if there are free parameters
					for(z in 1:length(setCollapseZero)){
						positionOfFirstCol<-setCollapseZero[z]
						startGridFirstPart<-data.frame(startGrid[,sequence(positionOfFirstCol - 1)])
						names(startGridFirstPart)<-names(startGrid)[sequence(positionOfFirstCol - 1)]
						if((1 + length(sequence(positionOfFirstCol - 1))) > length(names(startGrid))){
							startGridSecondPart<-c()
							startGrid<-cbind(startGridFirstPart,log(0))
						}else{
							startGridSecondPart<-startGrid[(1 + length(sequence(positionOfFirstCol - 1))):length(names(startGrid))]
							startGrid<-cbind(startGridFirstPart,log(0),startGridSecondPart)
						}
				
						names(startGrid)[positionOfFirstCol]<-paste("collapse_",setCollapseZero[z],sep="")
					}
				}
  			
				#If there are added events in addition to zero-fixed collapses, need to re-filter the grid to make sure
				#that added event times make sense in the context of the specified history
				if(!is.null(addedEventTime)){
					startGrid<-CreateStartGrid(startGrid,migrationIndividual=migrationArray[[modelID]],
						addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,
						startGridWithSetCollapseZero=TRUE)
					startGrid<-startGrid[[1]]
				}
			}

			#Stitch in a 1 to replace the first n0 (which is never free)	
			positionOfFirstN0 <- grep("n0multiplier", paramNames)[1]
			startGridFirstPart<-data.frame(startGrid[,sequence(positionOfFirstN0 - 1)])
			names(startGridFirstPart)<-names(startGrid)[sequence(positionOfFirstN0 - 1)]
			if((1 + length(sequence(positionOfFirstN0 - 1))) > length(names(startGrid))){
				startGridSecondPart<-c()
				startGrid<-cbind(startGridFirstPart,log(1))
			}else{
				if(length(grep("n0multiplier", paramNames)) <= 1){
					startGridSecondPart<-startGrid[(1 + length(sequence(positionOfFirstN0 - 1))):length(names(startGrid))]
				}else{
					startGridSecondPart<-startGrid[(2 + length(sequence(positionOfFirstN0 - 1))):length(names(startGrid))]
				}
				startGrid<-cbind(startGridFirstPart,log(1),startGridSecondPart)
			}
			names(startGrid)[positionOfFirstN0]<-"n0multiplier_1"
		}            
      
		if(is.null(nTreesGrid)){
			nTreesGrid<-10*nTrees #thinking here that want better estimate on the grid than in the heat of the search
		}
  	
		#Get and store AIC for each set of grid values (if checkpointing enabled, import and output AICs as necessary)
		if(!is.null(checkpointFile)){
			if(file.exists(checkpointFile)){
				if(as.numeric(strsplit(system(paste("ls -s ",checkpointFile,sep=""),intern=TRUE)," ")[[1]][1]) != 0){
					initial.AIC<-as.array(read.table(checkpointFile,stringsAsFactors=FALSE)[,1])
					startingPosition<-length(initial.AIC) + 1
				}else{
					initial.AIC<-c()
					startingPosition<-1
				}
			}else{
				initial.AIC<-c()
				startingPosition<-1
			}
			if(file.exists(checkpointFile)){
				unlink(checkpointFile)
			}
		}else{
			initial.AIC<-c()
			startingPosition<-1
		}

		if(!is.null(checkpointFile)){
			write.table(initial.AIC,file=checkpointFile,quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)
		}
	
		if(startingPosition <= nrow(startGrid)){
			for(t in startingPosition:nrow(startGrid)){
				currentAIC<-ReturnAIC(par=as.numeric(startGrid[t,]),migrationIndividual=migrationArray[[modelID]],
				badAIC=badAIC,maxParameterValue=maxParameterValue,nTrees=nTreesGrid,msPath=msPath,comparePath=comparePath,
				unresolvedTest=unresolvedTest,print.ms.string=print.ms.string,print.results=print.results,print.matches=print.matches,
				ncores=ncores,debug=debug,numReps=numReps,parameterBounds=parameterBounds,subsamplesPerGene=subsamplesPerGene,
				totalPopVector=totalPopVector,popAssignments=popAssignments,subsampleWeights.df=subsampleWeights.df,
				summaryFn=summaryFn,saveNoExtrap=saveNoExtrap,doSNPs=doSNPs,nEq=nEq,setCollapseZero=setCollapseZero,rm.n0=rm.n0,
				popScaling=popScaling,addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,
				optimization=optimization)
			
				initial.AIC<-append(initial.AIC,currentAIC)
				if(!is.null(checkpointFile)){
					write.table(currentAIC,file=checkpointFile,quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE,append=TRUE)
				}
			}
		}
		if(debug){
			if(rm.n0){
				print(cbind(AIC=initial.AIC,exp(startGrid[,-(grep("n0multiplier",colnames(startGrid)))])))
			}else{
				print(cbind(AIC=initial.AIC,exp(startGrid)))
			}
		}

		#Save the grid as object
		positionOfFirstN0 <- grep("n0multiplier", MsIndividualParameters(migrationArray[[modelID]]))[1]
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
    
   
		#Finish off with nloptr optimization   
		#Get starting parameters
		if(optimization == "nloptr"){
			prunedStartGrid1 <- as.matrix(startGrid[-grep("n0multiplier_1",names(startGrid))]) #Always toss the first n0multiplier as it is never free
			prunedStartGrid2 <- matrix(prunedStartGrid1[which(!is.na(thisGrid$AIC)),],ncol=ncol(prunedStartGrid1))
			colnames(prunedStartGrid2)<-colnames(prunedStartGrid1)
		
			for(rep in sequence(min(numReps, dim(startGrid)[1]))) {

				if(rep > nrow(prunedStartGrid2)){ #If no parameter estimates available from a gridSearch, take a random grid combination
					prunedStartGrid3<-unlist(prunedStartGrid1[order(initial.AIC)[sample(1:nrow(prunedStartGrid1),1)],])
				}else{ #else, pick the best parameter combination that hasn't already been used
					prunedStartGrid3<-unlist(prunedStartGrid2[rep,])
				}
	
			   #startingVals<-log(c(rlnorm(sum(grepl("collapse",paramNames)),1,1), rlnorm(-1+sum(grepl("n0multiplier",paramNames)),1,1), rbeta(sum(grepl("migration",paramNames)),shape1=1,shape2=3) )) #going to optimize in log space. Remove the first n0 parameter from optimization vector
			   startingVals <- prunedStartGrid3
	   
			   if(debug) {
				 print(startingVals) 
				}
	
				#Note, there are instances when nloptr is throwing an error (x0 returns NA objective). Need to troubleshoot this.
				searchResults<-nloptr(x0=startingVals, eval_f=ReturnAIC, opts=list("maxeval"=itnmax, "algorithm"="NLOPT_LN_SBPLX", "print_level"=1,
					maxtime=maxtime, maxeval=maxeval), migrationIndividual=migrationArray[[modelID]],
					badAIC=badAIC, maxParameterValue=maxParameterValue,
					nTrees=nTrees,msPath=msPath,comparePath=comparePath,unresolvedTest=unresolvedTest,ncores=ncores,
					print.ms.string=print.ms.string, print.results=print.results,print.matches=print.matches,
					debug=debug,numReps=numReps,parameterBounds=parameterBounds,subsamplesPerGene=subsamplesPerGene,summaryFn=summaryFn,
					totalPopVector=totalPopVector,subsampleWeights.df=subsampleWeights.df,popAssignments=popAssignments,
					saveNoExtrap=saveNoExtrap,doSNPs=doSNPs,nEq=nEq,setCollapseZero=setCollapseZero,rm.n0=rm.n0,popScaling=popScaling,
					addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,optimization=optimization)
				#If setCollapseZero is being used, convert optimized parameters to log(0)
				if(!is.null(setCollapseZero)){
					searchResults$solution[setCollapseZero]<-log(0)
				}
		
				#Keep best values
				if(!is.na(searchResults$objective)){
					if(searchResults$objective <= best.result.objective) {
						best.result<-searchResults
						best.result.objective<-searchResults$objective	
					}
				}else{ #If the AIC estimated = NA, construct dummy result
					best.result<-list("objective"=NA,"solution"=unname(startingVals))
				}
			}
		}
		best.resultGrid<-list(best.result,thisGrid)
		ifelse(return.all, return(best.resultGrid), return(best.result.objective))   
	}  
}


SearchContinuousModelSpaceRGenoud<-function(p, migrationArrayMap=NULL, migrationArray, popAssignments, badAIC=100000000000000, 
	maxParameterValue=20, nTrees=2e5 ,msPath="ms",comparePath=system.file("extdata", "comparecladespipe.pl", package="phrapl"),
	subsampleWeights.df=NULL,
	unresolvedTest=TRUE,print.ms.string=FALSE, ncores=1,print.results=TRUE,print.matches=FALSE,debug=FALSE,method="nlminb",
	itnmax=NULL, return.all=FALSE, maxtime=0, maxeval=0, parameterBounds=list(minCollapseTime=0.1,
	minCollapseRatio=0,minN0Ratio=0.1,minGrowth=0.1,minGrowthRatio=0.1,minMigrationRate=0.05,minMigrationRatio=0.1), 
	numReps=0, startGrid=startGrid, collapseStarts=c(0.30,0.58,1.11,2.12,4.07,7.81,15.00), n0Starts=c(0.1,0.5,1,2,4), 
	growthStarts=c(0.30,0.58,1.11,2.12,4.07,7.81,15.00),migrationStarts=c(0.10,0.22,0.46,1.00,2.15,4.64), gridSave=NULL,
	gridSaveFile=NULL,subsamplesPerGene=1,totalPopVector,summaryFn="mean",saveNoExtrap=FALSE,doSNPs=FALSE,nEq=100,
	setCollapseZero=NULL,rm.n0=TRUE,popScaling=NULL,checkpointFile=NULL,addedEventTime=NULL,addedEventTimeAsScalar=TRUE,
	genoudPopSize=25,numGridStartVals=25,solutionTolerance=1,skipGrid=FALSE, ...) {

	if(!is.null(migrationArrayMap)){
  		modelID<-ReturnModel(p,migrationArrayMap)
  	}else{
  		modelID<-p
  	}
  	best.result <- c()
  	best.result.objective <- badAIC
  	if(is.na(modelID)) {
		return(badAIC)
  	}else{
		paramNames<-MsIndividualParameters(migrationArray[[modelID]])
		if(is.null(startGrid)) { #need to make our own grid
			startingVectorList<-list()
			nameCount<-1
		
			#Calculate number of collapse parameters to estimate (subtracting off those that are set to zero)
			numCollapse <- sum(grepl("collapse",paramNames)) - length(setCollapseZero)
			if(numCollapse < 0){ #If for some reason there are more fixed parameters than possible, reduce them
				setCollapseZero<-sequence(sum(grepl("collapse",paramNames)))
				numCollapse <- sum(grepl("collapse",paramNames)) - length(setCollapseZero)
			}
			if(!is.null(setCollapseZero)){ #If setCollapseZero is nonzero
				collapses2estimate<-grep("collapse",paramNames,value=TRUE)[which(!grep("collapse",paramNames) %in% setCollapseZero)]
			}else{
				collapses2estimate<-grep("collapse",paramNames,value=TRUE)
			}
			if(numCollapse > 0) {
				for (i in sequence(numCollapse)) {
					startingVectorList<-append(startingVectorList, list(collapseStarts))
					names(startingVectorList)[nameCount]<-collapses2estimate[i]
					nameCount<-nameCount + 1
				} 
			}
		
			#Calculate number of n0 parameters (if just one, subtract it off)
			if(sum(grepl("n0multiplier",paramNames)) == 1){
			   numn0 <- sum(grepl("n0multiplier",paramNames)) - 1 # -1 since one is fixed to be 1
			}else{
			   numn0 <- sum(grepl("n0multiplier",paramNames))
			}
			if (numn0 > 0) {
				for (i in sequence(numn0)) {
					startingVectorList<-append(startingVectorList, list(n0Starts))
					names(startingVectorList)[nameCount]<-grep("n0multiplier",paramNames,value=TRUE)[i]
					nameCount<-nameCount + 1
				} 
			}	
		
			#Calculate number of growth parameters
			numGrowth <- sum(grepl("growth",paramNames))
			if (numGrowth > 0) {
				for (i in sequence(numGrowth)) {
					startingVectorList<-append(startingVectorList, list(growthStarts))
					names(startingVectorList)[nameCount]<-grep("growth",paramNames,value=TRUE)[i]
					nameCount<-nameCount + 1
				} 
			}
		
			#Calculate number of migration parameters
			numMigration <- sum(grepl("migration",paramNames))
			if (numMigration > 0) {
				for (i in sequence(numMigration)) {
					startingVectorList<-append(startingVectorList, list(migrationStarts))
					names(startingVectorList)[nameCount]<-grep("migration",paramNames,value=TRUE)[i]
					nameCount<-nameCount + 1
				} 
			}
		
			#If some collapses are set to zero, there are collapses in the model, and also added events,
			#wait to filter the grid for added events until after zero-fixed collapses are stitched in below
			if(!is.null(setCollapseZero) && sum(grepl("collapse",paramNames)) != 0 && !is.null(addedEventTime)){
				startGrid<-CreateStartGrid(fineGrid=startingVectorList,migrationIndividual=migrationArray[[modelID]])
				startGrid<-startGrid[[1]] #default grid shouldn't be list (as there is always one grid)
		
			#Otherwise, filter for any added events now
			}else{
				startGrid<-CreateStartGrid(fineGrid=startingVectorList,migrationIndividual=migrationArray[[modelID]],
					addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar)
				startGrid<-startGrid[[1]] #default grid shouldn't be list (as there is always one grid)
			}
		
			#If some collapses are set to zero (and there are collapses in the model), stitch these values into the startGrid
			if(!is.null(setCollapseZero) && sum(grepl("collapse",paramNames)) != 0){
			
				#If there are only fixed parameters...
				if(length(startGrid) == 0){
					startGrid<-data.frame(matrix(rep(log(0),length(setCollapseZero)),nrow=1,ncol=length(setCollapseZero)))
					for(g in 1:length(setCollapseZero)){
						names(startGrid)[g]<-paste("collapse_",setCollapseZero[g],sep="")
					}
				}else{
			
					#Otherwise, if there are free parameters
					for(z in 1:length(setCollapseZero)){
						positionOfFirstCol<-setCollapseZero[z]
						startGridFirstPart<-data.frame(startGrid[,sequence(positionOfFirstCol - 1)])
						names(startGridFirstPart)<-names(startGrid)[sequence(positionOfFirstCol - 1)]
						if((1 + length(sequence(positionOfFirstCol - 1))) > length(names(startGrid))){
							startGridSecondPart<-c()
							startGrid<-cbind(startGridFirstPart,log(0))
						}else{
							startGridSecondPart<-startGrid[(1 + length(sequence(positionOfFirstCol - 1))):length(names(startGrid))]
							startGrid<-cbind(startGridFirstPart,log(0),startGridSecondPart)
						}
				
						names(startGrid)[positionOfFirstCol]<-paste("collapse_",setCollapseZero[z],sep="")
					}
				}
			
				#If there are added events in addition to zero-fixed collapses, need to re-filter the grid to make sure
				#that added event times make sense in the context of the specified history
				if(!is.null(addedEventTime)){
					startGrid<-CreateStartGrid(startGrid,migrationIndividual=migrationArray[[modelID]],
						addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,
						startGridWithSetCollapseZero=TRUE)
					startGrid<-startGrid[[1]]
				}
			}

			#Stitch in a 1 to replace the first n0 (which is never free)	
			positionOfFirstN0 <- grep("n0multiplier", paramNames)[1]
			startGridFirstPart<-data.frame(startGrid[,sequence(positionOfFirstN0 - 1)])
			names(startGridFirstPart)<-names(startGrid)[sequence(positionOfFirstN0 - 1)]
			if((1 + length(sequence(positionOfFirstN0 - 1))) > length(names(startGrid))){
				startGridSecondPart<-c()
				startGrid<-cbind(startGridFirstPart,log(1))
			}else{
				if(length(grep("n0multiplier", paramNames)) <= 1){
					startGridSecondPart<-startGrid[(1 + length(sequence(positionOfFirstN0 - 1))):length(names(startGrid))]
				}else{
					startGridSecondPart<-startGrid[(2 + length(sequence(positionOfFirstN0 - 1))):length(names(startGrid))]
				}
				startGrid<-cbind(startGridFirstPart,log(1),startGridSecondPart)
			}
			names(startGrid)[positionOfFirstN0]<-"n0multiplier_1"
		}
		
	
		#Get and store AIC for each set of grid values (if checkpointing enabled, import and output AICs as necessary)
		if(!is.null(checkpointFile)){
			if(file.exists(checkpointFile)){
				if(as.numeric(strsplit(system(paste("ls -s ",checkpointFile,sep=""),intern=TRUE)," ")[[1]][1]) != 0){
					initial.AIC<-as.array(read.table(checkpointFile,stringsAsFactors=FALSE)[,1])
					startingPosition<-length(initial.AIC) + 1
				}else{
					initial.AIC<-c()
					startingPosition<-1
				}
			}else{
				initial.AIC<-c()
				startingPosition<-1
			}
			if(file.exists(checkpointFile)){
				unlink(checkpointFile)
			}
		}else{
			initial.AIC<-c()
			startingPosition<-1
		}

		if(!is.null(checkpointFile)){
			write.table(initial.AIC,file=checkpointFile,quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)
		}
	
		#If you want to skip the grid, and just move on to optimization
		
		if(startingPosition <= nrow(startGrid)){
			for(t in startingPosition:nrow(startGrid)){
				if(skipGrid != TRUE){ #if for optimization you want to skip the grid, just fill it with NAs
					currentAIC<-ReturnAIC(par=as.numeric(startGrid[t,]),migrationIndividual=migrationArray[[modelID]],
					badAIC=badAIC,maxParameterValue=maxParameterValue,nTrees=nTrees,msPath=msPath,comparePath=comparePath,
					unresolvedTest=unresolvedTest,print.ms.string=print.ms.string,print.results=print.results,print.matches=print.matches,
					ncores=ncores,debug=debug,numReps=numReps,parameterBounds=parameterBounds,subsamplesPerGene=subsamplesPerGene,
					totalPopVector=totalPopVector,popAssignments=popAssignments,subsampleWeights.df=subsampleWeights.df,
					summaryFn=summaryFn,saveNoExtrap=saveNoExtrap,doSNPs=doSNPs,nEq=nEq,setCollapseZero=setCollapseZero,rm.n0=rm.n0,
					popScaling=popScaling,addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,
					optimization=optimization)
	
					initial.AIC<-append(initial.AIC,currentAIC)
				}else{
					initial.AIC<-append(initial.AIC,NA)
				}
				
				if(!is.null(checkpointFile)){
					write.table(currentAIC,file=checkpointFile,quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE,append=TRUE)
				}
			}
		}		
		
		if(debug){
			if(rm.n0){
				print(cbind(AIC=initial.AIC,exp(startGrid[,-(grep("n0multiplier",colnames(startGrid)))])))
			}else{
				print(cbind(AIC=initial.AIC,exp(startGrid)))
			}
		}

		#Save the grid as object
		positionOfFirstN0 <- grep("n0multiplier", MsIndividualParameters(migrationArray[[modelID]]))[1]
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

		#Get starting parameters for genoud, starting with the best parameter values from a gridSearch
		prunedStartGrid1 <- data.frame(as.matrix(startGrid[-grep("n0multiplier_1",names(startGrid))]))
		prunedStartGrid2 <- prunedStartGrid1[which(!is.na(thisGrid$AIC)),]
		if(nrow(prunedStartGrid2) == 0){ #If no parameter estimates available from a gridSearch, take the grid median
			prunedStartGrid2<-matrix(sapply(exp(as.data.frame(prunedStartGrid1)),median),nrow=1,ncol=ncol(prunedStartGrid2))
			colnames(prunedStartGrid2)<-colnames(prunedStartGrid1)
			prunedStartGrid2<-log(prunedStartGrid2)
		}else{
			if(nrow(prunedStartGrid2) >= numGridStartVals){ #Take the best values from grid search, if available
				prunedStartGrid2<-prunedStartGrid2[1:numGridStartVals,]
			}else{ #Or take all that's available
				prunedStartGrid2<-prunedStartGrid2[1:nrow(prunedStartGrid2),]
			}
		}
		
		#Add in other random grid points if there aren't enough values with estimable AICs
		#if skipGrid = TRUE, then all starting values will be randomly selected
		prunedStartGridNAs<-prunedStartGrid1[which(is.na(thisGrid$AIC)),]
		if(nrow(prunedStartGrid2) < numGridStartVals){
			if(nrow(prunedStartGridNAs) >= numGridStartVals){
				randomValues<-sort(sample(1:nrow(prunedStartGridNAs),(numGridStartVals - nrow(prunedStartGrid2))))
				prunedStartGrid2<-rbind(prunedStartGrid2,prunedStartGridNAs[randomValues,])
			}else{
				prunedStartGrid2<-rbind(prunedStartGrid2,prunedStartGridNAs)
			}
		}
		
		#Do genoud search
		nvars <- dim(prunedStartGrid2)[2]
		prunedStartGrid2<-as.matrix(prunedStartGrid2)
		searchResults<-genoud(fn=ReturnAIC,nvars= nvars,pop.size=genoudPopSize,starting.values=prunedStartGrid2, 
			Domains=matrix(c(rep(log(0.01), nvars), rep(log(maxParameterValue), nvars)), ncol=2, byrow=FALSE),boundary.enforcement=1, 
			gradient.check=FALSE,hessian=FALSE,solution.tolerance=solutionTolerance,wait.generations=5,migrationIndividual=migrationArray[[modelID]], 
			badAIC=badAIC, maxParameterValue=maxParameterValue,
			nTrees=nTrees,msPath=msPath,comparePath=comparePath,unresolvedTest=unresolvedTest,ncores=ncores,
			print.ms.string=print.ms.string, print.results=print.results,print.matches=print.matches,
			debug=debug,numReps=numReps,parameterBounds=parameterBounds,subsamplesPerGene=subsamplesPerGene,summaryFn=summaryFn,
			totalPopVector=totalPopVector,subsampleWeights.df=subsampleWeights.df,popAssignments=popAssignments,
			saveNoExtrap=saveNoExtrap,doSNPs=doSNPs,nEq=nEq,setCollapseZero=setCollapseZero,rm.n0=rm.n0,popScaling=popScaling,
			addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,optimization=optimization)

			#If setCollapseZero is being used, convert optimized parameters to log(0)
			if(!is.null(setCollapseZero)){
				searchResults$solution[setCollapseZero]<-log(0)
			}
		
			#Keep best values and grid
			best.resultGrid<-list(searchResults,thisGrid)
	
		return(best.resultGrid)    
	}
}

#This function was formally known as "ExhaustiveSearchNLoptr", and still has the capabilities of running a heuristic
#nloptr search. However, since we are currently focusing on a grid search, we have changed the name of the function
GridSearch<-function(modelRange=c(1:length(migrationArray)),migrationArrayMap=NULL,migrationArray,popAssignments, 
		badAIC=100000000000000,maxParameterValue=20,nTrees=1e5 ,msPath="ms",comparePath=system.file("extdata", "comparecladespipe.pl", package="phrapl"),
		observedTrees,subsampleWeights.df=NULL, doWeights=TRUE, unresolvedTest=TRUE,
		print.ms.string=FALSE,print.results=TRUE,print.matches=FALSE,debug=FALSE,method="nlminb",itnmax=NULL,
		ncores=1,results.file=NULL,maxtime=0, maxeval=0,return.all=TRUE, numReps=0,startGrid=NULL,
		collapseStarts=c(0.30,0.58,1.11,2.12,4.07,7.81,15.00), n0Starts=c(0.1,0.5,1,2,4), 
		growthStarts=c(0.30,0.58,1.11,2.12,4.07,7.81,15.00),migrationStarts=c(0.1,0.22,0.46,1,2.15,4.64),
		subsamplesPerGene=1,totalPopVector=NULL,summaryFn="mean",saveNoExtrap=FALSE,doSNPs=FALSE,
		nEq=100,setCollapseZero=NULL,dAIC.cutoff=2,rm.n0=TRUE,popScaling=NULL,checkpointFile=NULL,
		parameter.ambiguous=FALSE,addedEventTime=NULL,addedEventTimeAsScalar=TRUE,optimization="grid",
		genoudPopSize=25,numGridStartVals=25,solutionTolerance=1,skipGrid=FALSE,...){		
	#If multiple n0multiplier values exists in the migrationArray, then make rm.n0 FALSE, else, leave it as specified
	n0Values<-unlist(migrationArray)
	if(length(unique(n0Values[grep("n0multiplier",names(n0Values))][!is.na(n0Values[grep("n0multiplier",names(n0Values))])])) > 1){
		rm.n0<-FALSE
	}	
	
	#If no popScaling defined, assume same scalar across loci
	if(is.null(popScaling)) {
		popScaling <- rep(1, length(observedTrees[[1]]))
	}else{
		#If popScaling defined, expand popScaling to repeat across subsamples
		popScaling<-rep(popScaling,subsamplesPerGene)
	}
	
	if(length(modelRange) != length(migrationArray)) { #need to look at only particular rows
		migrationArray<-migrationArray[modelRange]
		if(!is.null(migrationArrayMap)){
			migrationArrayMap<-GenerateMigrationArrayMapTrunc(migrationArrayMap,modelRange)
		}
	}
	if(is.null(subsampleWeights.df) && doWeights) {
   	subsampleWeights.df <- GetPermutationWeightsAcrossSubsamples(popAssignments, observedTrees)
	}
  	#Prepare temporary tree and assign files
	for(k in 1:length(popAssignments)){
		write.table(CreateAssignment.df(popAssignments[[k]]),file=paste(tempdir(),"/assign",k,".txt",sep=""),
			quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,append=FALSE)
		for(scaling.index in sequence(length(unique(popScaling)))) {
			if(file.exists(paste(tempdir(),"/observed",k,".scaling.", popScaling[scaling.index],".tre",sep=""))){
				unlink(paste(tempdir(),"/observed",k,".scaling.", popScaling[scaling.index],".tre",sep=""))
			}
		}
	}
	for(f in 1:length(observedTrees)){ #for each popAssignment
		if (class(observedTrees[[f]]) != "multiPhylo"){
        	observedTrees[[f]] <- c(observedTrees[[f]])
    	}
		for(m in 1:length(observedTrees[[f]])){ #for each tree
			write.tree(observedTrees[[f]][[m]],file=paste(tempdir(),"/observed",f,".scaling.", popScaling[m], ".tre",sep=""),append=TRUE)
		}		
	}
 
  AIC.values<-rep(NA,length(migrationArray))
  gridList<-list() #for storing model grids
  results.list<-list()
  for (i in sequence(length(migrationArray))){
  	if(!is.null(migrationArrayMap)){
  		p<-c(migrationArrayMap$collapseMatrix.number[i], migrationArrayMap$n0multiplierMap.number[i], 
  			migrationArrayMap$growthMap.number[i], migrationArrayMap$migrationArray.number[i])
  	}else{
  		p<-i
  	}
  	if(return.all) {
 		if(!is.null(startGrid)){
 			currentStartGrid<-startGrid[[i]] #imported startGrid must be a list
 		}else{
 			currentStartGrid=NULL
 		}

  		result.indiv<-NULL
  		if(optimization != "genoud"){
			try(result.indiv<-SearchContinuousModelSpaceNLoptr(p,migrationArrayMap,migrationArray,popAssignments,badAIC=badAIC,
				maxParameterValue=maxParameterValue,nTrees=nTrees,msPath=msPath,comparePath=comparePath,unresolvedTest=unresolvedTest,
				print.ms.string=print.ms.string,print.results=print.results,print.matches=print.matches,debug=debug,
				method=method,itnmax=itnmax, maxtime=maxtime,ncores=ncores,maxeval=maxeval,return.all=return.all,numReps=numReps,
				startGrid=currentStartGrid,gridSave=NULL,collapseStarts=collapseStarts,n0Starts=n0Starts,growthStarts=growthStarts,
				migrationStarts=migrationStarts,subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,
				subsampleWeights.df=subsampleWeights.df,summaryFn=summaryFn,saveNoExtrap=saveNoExtrap,doSNPs=doSNPs,nEq=nEq,
				setCollapseZero=setCollapseZero,rm.n0=rm.n0,popScaling=popScaling,checkpointFile=checkpointFile,
				addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,optimization=optimization, ...))
		}else{
			try(result.indiv<-SearchContinuousModelSpaceRGenoud(p,migrationArrayMap,migrationArray,popAssignments,badAIC=badAIC,
				maxParameterValue=maxParameterValue,nTrees=nTrees,msPath=msPath,comparePath=comparePath,unresolvedTest=unresolvedTest,
				print.ms.string=print.ms.string,print.results=print.results,print.matches=print.matches,debug=debug,
				method=method,itnmax=itnmax, maxtime=maxtime,ncores=ncores,maxeval=maxeval,return.all=return.all,numReps=numReps,
				startGrid=currentStartGrid,gridSave=NULL,collapseStarts=collapseStarts,n0Starts=n0Starts,growthStarts=growthStarts,
				migrationStarts=migrationStarts,subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,
				subsampleWeights.df=subsampleWeights.df,summaryFn=summaryFn,saveNoExtrap=saveNoExtrap,doSNPs=doSNPs,nEq=nEq,
				setCollapseZero=setCollapseZero,rm.n0=rm.n0,popScaling=popScaling,checkpointFile=checkpointFile,
				addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,genoudPopSize=genoudPopSize,
				numGridStartVals=numGridStartVals,solutionTolerance=solutionTolerance,optimization=optimization,skipGrid=skipGrid, ...))
		}
	
  		if(rm.n0){
  			result.indiv[[2]]<-result.indiv[[2]][,-(grep("n0multiplier",colnames(result.indiv[[2]])))]
  		}
  		
		gridList[[length(gridList)+1]]<-result.indiv[[2]] #make list of model grids

  		if(!is.null(result.indiv[[1]])) {
  			if(optimization == "nloptr"){
  				AIC.values[i]<-result.indiv[[1]]$objective	
  			}
  			if(optimization == "genoud"){
  				AIC.values[i]<-result.indiv[[1]]$value	
  			}
  		}
  		results.list<-append(results.list, list(result.indiv[[1]]))

  	}else{
  		if(optimization != "genoud"){
			try(AIC.values[i]<-SearchContinuousModelSpaceNLoptr(p,migrationArrayMap,migrationArray,popAssignments,badAIC=badAIC,
			maxParameterValue=maxParameterValue, nTrees=nTrees, msPath=msPath,comparePath=comparePath,ncores=ncores,
			unresolvedTest=unresolvedTest,print.ms.string=print.ms.string,print.results=print.results,print.matches=print.matches,
			debug=debug,method=method,itnmax=itnmax, maxtime=maxtime, maxeval=maxeval, return.all=return.all,numReps=numReps,
			startGrid=currentStartGrid,collapseStarts=collapseStarts,n0Starts=n0Starts,growthStarts=growthStarts, migrationStarts=migrationStarts,
			gridSave=NULL,subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,subsampleWeights.df=subsampleWeights.df,
			summaryFn=summaryFn,saveNoExtrap=saveNoExtrap,doSNPs=doSNPs,nEq=nEq,setCollapseZero=setCollapseZero,rm.n0=rm.n0,popScaling=popScaling,
			checkpointFile=checkpointFile,addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,optimization=optimization,...))
		}else{
			try(AIC.values[i]<-SearchContinuousModelSpaceRGenoud(p,migrationArrayMap,migrationArray,popAssignments,badAIC=badAIC,
			maxParameterValue=maxParameterValue, nTrees=nTrees, msPath=msPath,comparePath=comparePath,ncores=ncores,
			unresolvedTest=unresolvedTest,print.ms.string=print.ms.string,print.results=print.results,print.matches=print.matches,
			debug=debug,method=method,itnmax=itnmax, maxtime=maxtime, maxeval=maxeval, return.all=return.all,numReps=numReps,
			startGrid=currentStartGrid,collapseStarts=collapseStarts,n0Starts=n0Starts,growthStarts=growthStarts, migrationStarts=migrationStarts,
			gridSave=NULL,subsamplesPerGene=subsamplesPerGene,totalPopVector=totalPopVector,subsampleWeights.df=subsampleWeights.df,
			summaryFn=summaryFn,saveNoExtrap=saveNoExtrap,doSNPs=doSNPs,nEq=nEq,setCollapseZero=setCollapseZero,rm.n0=rm.n0,popScaling=popScaling,
			checkpointFile=checkpointFile,addedEventTime=addedEventTime,addedEventTimeAsScalar=addedEventTimeAsScalar,
			genoudPopSize=genoudPopSize,numGridStartVals=numGridStartVals,solutionTolerance=solutionTolerance,optimization=optimization,
			skipGrid=skipGrid, ...))
		}
  	}

#  	 	print(c(i, length(migrationArray), i/length(migrationArray), AIC.values[i]))

  		if(!is.null(results.file)) {
			save(list=ls(), file=results.file)
		}
	}

  	#Toss temporary tree and assign files
# 	for(k in 1:length(popAssignments)){
# 		unlink(c(paste(tempdir(),"/assign",k,".txt",sep=""),paste(tempdir(),"/observed",k,".tre",sep=""),paste(tempdir(), "/mstrees.txt", sep="")))
# 	}
#	print(results.list)


	#Save table of best models
	if(optimization == "grid"){
		overall.results<-ExtractGridAICs(result=gridList,migrationArray=migrationArray,setCollapseZero=setCollapseZero)
	}else{
		overall.results<-ExtractOptimizationAICs(result=results.list,migrationArray=migrationArray,setCollapseZero=setCollapseZero,
			optimization=optimization)
	}
	overall.results$models <- modelRange #so if we've done subsampling, use right model numbers

	#If the grid list is filled with NAs (i.e., AIC couldn't be estimated given that all parameter
	#combinations resulted in too many trees with zero matches), for all models, give an error
	if(skipGrid != TRUE){
		if(length(which(!is.na(overall.results$AIC))) == 0){
			warning("Error: There are not enough trees being simulated to estimate a log-likelihood for any of the models. Increase nTrees or reduce the number of tips in the subsampled trees\n")
			return()
		}
	}
	
	#If the grid list for only one of several analyzed models is filled with NAs, just toss results for those models.
	if(length(which(is.na(overall.results$AIC))) > 0){
		whichModelsToKeep<-which(!is.na(overall.results$AIC))
		modelIDTossed<-overall.results[-whichModelsToKeep,]$models
		if(length(modelIDTossed) > 1){
			modelIDTossed<-paste(modelIDTossed,collapse=", ")
		}
		overall.results<-overall.results[whichModelsToKeep,]
		gridList<-gridList[whichModelsToKeep]
		
		warning(paste("There are not enough trees being simulated to estimate a log-likelihood for model(s) ",
			modelIDTossed,". Increase nTrees or reduce the number of tips in the subsampled trees\n"))
	}

		
	####Get parameters using the old ambiguous method (ExtractGridParameters)
	if(parameter.ambiguous==TRUE){	
		#Save parameter estimates and parameter indexes to tables
		if(numReps==0){
			parameters<-ExtractGridParameters(migrationArray=migrationArray,result=gridList,
				popVector=popAssignments[[1]],dAIC.cutoff=dAIC.cutoff)
		}else{
			parameters<-ExtractParameters(migrationArray=migrationArray,result=results.list,
				popVector=popAssignments[[1]])
		}
	
		#Add all results to list
		if(numReps == 0 && optimization == "nloptr"){
			results.final<-list("AIC.Grid"=gridList,"overall.results"=overall.results,
				"parameters"=parameters[[1]],"parameterIndexes"=parameters[[2]])
		}else{
			results.final<-list("search.results"=results.list,"AIC.Grid"=gridList,
				"overall.results"=overall.results,"parameters"=parameters[[1]],"parameterIndexes"=parameters[[2]])
		}
	}else{
		
		#Get parameters using new unambiguous method (ExtractUnambiguousGridParameters)		
		if(optimization == "grid"){
			##Get unambiguous parameters (grid)
			parameters<-ExtractUnambiguousGridParameters(overall.results=overall.results,gridList=gridList,
				migrationArray=migrationArray,sortParameters=TRUE,sortModelsAIC=TRUE)
			parameters<-parameters[order(parameters$models),]
			parametersOnly<-data.frame(matrix(as.matrix(parameters[,-c(1:2)]),nrow=nrow(parameters)))
			colnames(parametersOnly)<-colnames(parameters)[-c(1:2)]
			overall.results<-cbind(overall.results,parametersOnly)
			
			##Concatenate overall.results and parameters
			results.final<-list("AIC.Grid"=gridList,"overall.results"=overall.results)
		}else{
			##Get unambiguous parameters (optimization)
			parameters<-ExtractUnambiguousOptimizationParameters(overall.results=overall.results,
				results.list=results.list,migrationArray=migrationArray,sortParameters=TRUE,
				sortModelsAIC=TRUE,optimization=optimization)
			parameters<-parameters[order(parameters$models),]
			parametersOnly<-data.frame(matrix(as.matrix(parameters[,-c(1:2)]),nrow=nrow(parameters)))
			colnames(parametersOnly)<-colnames(parameters)[-c(1:2)]
			overall.results<-cbind(overall.results,parametersOnly)
			
			##Concatenate overall.results and parameters
			results.final<-list("search.results"=results.list,"AIC.Grid"=gridList,
				"overall.results"=overall.results)
		}
	}

	ifelse(return.all, return(results.final), return(AIC.values))     
}