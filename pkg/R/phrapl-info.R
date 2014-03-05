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
	return(TRUE)
}

#maxParameterValue prevents MS from going nuts (not finishing) with really high migration or other rates
ReturnAIC<-function(par,popVector,migrationIndividual,nTrees=1,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, print.results=FALSE, print.ms.string=FALSE, debug=FALSE, badAIC=100000000000000, maxParameterValue=100, parameterBounds=list(minCollapseTime=0.1, minCollapseRatio=0, minN0Ratio=0.1, minMigrationRate=0.05, minMigrationRatio=0.1)) {
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
  likelihoodVector<-PipeMS(popVector=popVector,migrationIndividual=migrationIndividual,parameterVector=parameterVector,nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest,print.ms.string=print.ms.string, debug=debug)
  lnLValue<-ConvertOutputVectorToLikelihood(likelihoodVector, nTrees=nTrees, probOfMissing=1/howmanytrees(sum(popVector)))
  AICValue<-2*(-lnLValue + KAll(migrationIndividual))
  if(print.results) {
    resultsVector<-c(AICValue,lnLValue,parameterVector)
    names(resultsVector)<-c("AIC","lnL",MsIndividualParameters(migrationIndividual))
    print(resultsVector)
  }
  return(AICValue)
}
