library(partitions)

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

setMaxK<-function(popVector,nTrees=1,taxaPerParam=5) {
  maxK<-max(1,floor(sum(popVector)/(nTrees*taxaPerParam))) 
  return(maxK)
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
    for (i in 1: max(collapseMatrix,na.rm=TRUE)) {
      parameterList<-append(parameterList,paste("collapse_",i,sep="")) 
    }
  }
  for (i in 1:max(n0multiplierMap,na.rm=TRUE)) {
    parameterList<-append(parameterList,paste("n0multiplier_",i,sep="")) 
  }
  for (i in 1:max(migrationArray,na.rm=TRUE)) {
    parameterList<-append(parameterList,paste("migration_",i,sep="")) 
  }
  return(parameterList)
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
     initialN0multipler<-paste(initialN0multipler,"-n",i,n0multiplierParameters[n0multiplierMap[i,1] ],sep=" ")
    }
 
    initialMigration<-"-ma"
    for (i in 1:nPop) {
      for (j in 1:nPop) {
       if (i==j) {
          initialMigration<-paste(initialMigration, "x", sep=" ") 
       }
       else {
         initialMigration<-paste(initialMigration, migrationParameters[migrationArray[i,j,1] ],sep=" ")
       }
      }
    }
    msString<-paste(msString,initialN0multipler,initialMigration,sep=" ")
    
    #is there a collapse in the first gen?
    if (max(collapseMatrix[,1],na.rm=TRUE)>0) { #remember that everything goes to the lowest index population with state ==1
      initialCollapse<-""
      popsToCollapse<-which(collapseMatrix[,1]==1)
      for (i in 2:length(popsToCollapse)) {
        initialCollapse<-paste(initialCollapse, "-ej", collapseParameters[1], popsToCollapse[i], popsToCollapse[1]  ,sep=" ")       
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
          msString<-paste(msString,"-en",collapseTime,n0multiplierPos,n0multiplierParameters[n0multiplierMap[n0multiplierPos,generation] ],sep=" ")
        }

        for (fromPos in 1:nPop) {
          for (toPos in 1:nPop) {
            migrationArrayValue<-migrationArray[fromPos,toPos,generation]
           if (!is.na(migrationArrayValue) ) {
            msString<-paste(msString,"-em",collapseTime,fromPos,toPos,migrationParameters[migrationArrayValue],sep=" ")              
           }
          }
        }
        
        
       #is there a collapse in this gen?
       if (max(collapseMatrix[,generation],na.rm=TRUE)>0) { #remember that everything goes to the lowest index population with state ==1
        initialCollapse<-""
        popsToCollapse<-which(collapseMatrix[,generation]==1)
        for (i in 2:length(popsToCollapse)) {
          initialCollapse<-paste(initialCollapse, "-ej", collapseParameters[generation], popsToCollapse[i], popsToCollapse[1]  ,sep=" ")       
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
			#all those in the lastGen in state 1 were merged into one population, id'ed by the one with the lowest id
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
				rawIntervals<-rawIntervals[,colMax(rawIntervals)<=1] #we're okay with having all zeros: no population collapse
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
	firstIntervals<-firstIntervals[,colMax(firstIntervals)<=1] #we're okay with having all zeros: no population collapse
	firstIntervals<-firstIntervals[,colCountIf(firstIntervals,1)!=1] #we want intervals that have all zeros or at least two ones. What does it mean to have one lineage coalesce with itself?
	popIntervalsList<-list()
	for (i in 1:dim(firstIntervals)[2]) {
		popIntervalsList[[i]]<-popinterval(as.matrix(firstIntervals[,i],ncol=1))
	}
	popIntervalsList<-completeIntervals(updateCompletes(popIntervalsList))
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
generateN0multiplierIndividuals<-function(popVector,popIntervalsList=generateIntervals(popVector,maxK),maxK) {
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
					thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we've used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
				}
				n0multiplierIndividualsList[[length(n0multiplierIndividualsList)+1]]<-n0multiplierindividual(popIntervalsList[[i]]$collapseMatrix, popIntervalsList[[i]]$complete, n0multiplierMap)
			}
		}
	}
	return(n0multiplierIndividualsList)
}

kN0multiplierMap<-function(n0multiplierMap) {
  return(max(n0multiplierMap,na.rm=TRUE)) 
}

#now we will generate all possible assignments of pairwise migration. Again, we want to keep the total number of free parameters (times, n0multipliers, migration rates) under our chosen max
#allow a model where migrations change anywhere along branch, or only at coalescent nodes? The problem with the latter is that you can't fit some reasonable models: i.e., two populations persisting through time. Problem with the former is parameter space
generateMigrationIndividuals<-function(popVector,n0multiplierIndividualsList=generateN0multiplierIndividuals(popVector,popIntervalsList=generateIntervals(popVector,maxK),maxK), maxK, verbose=FALSE) {
	migrationIndividualsList<-list()
	for (i in 1:length(n0multiplierIndividualsList)) {
    if (verbose==TRUE) {
      print(paste("doing ",i,"/",length(n0multiplierIndividualsList)))
    }
		collapseMatrix<-n0multiplierIndividualsList[[i]]$collapseMatrix
    n0multiplierMap<-n0multiplierIndividualsList[[i]]$n0multiplierMap
		numFinalPops<-dim(collapseMatrix)[1]
    numSteps<-dim(collapseMatrix)[2]
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
		possibleMappings<-compositions(numPairs)

   for (mappingIndex in 1:dim(possibleMappings)[2]) {
			thisMapping<-possibleMappings[,mappingIndex]
			if ((length(which(thisMapping>0)) + kCollapseMatrix(collapseMatrix) + kN0multiplierMap(n0multiplierMap) )<=maxK) { #only do it on those mappings that have not too many free parameters
				migrationArray<-migrationTemplate
				whichPositions <- which(migrationArray==1)
				for (positionIndex in 1:length(whichPositions)) {
					position=whichPositions[positionIndex]
					paramPosition<-which(thisMapping>0)[1]
					migrationArray[position]=paramPosition #the position of the first parameter
					thisMapping[paramPosition]=thisMapping[paramPosition]-1 #now we've used up one of those parameters. If there are no more intervals assigned that parameter, it will drop to 0
				}
				migrationIndividualsList[[length(migrationIndividualsList)+1]]<-migrationindividual(n0multiplierIndividualsList[[i]]$collapseMatrix, n0multiplierIndividualsList[[i]]$complete, n0multiplierMap, migrationArray)
			}
		}		
	}
  return(migrationIndividualsList)
}

loadMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation="/usr/local/bin/ms") {
  msCallInfo<-createMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
  geneTrees<-system(paste(msLocation,msCallInfo$nsam,msCallInfo$nreps,msCallInfo$opts," | grep ';'"),intern=TRUE)
  return(geneTrees)
}

saveMS<-function(popVector,migrationIndividual,parameterVector,nTrees=1,msLocation="/usr/local/bin/ms",file="sim.tre") {
  msCallInfo<-createMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
  returnCode<-system(paste(msLocation,msCallInfo$nsam,msCallInfo$nreps,msCallInfo$opts," | grep ';' >",file),intern=FALSE)
  return(returnCode)
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

convertOutputVectorToLikelihood<-function(outputVector,nTrees) {
  outputVector<-as.numeric(outputVector)
  outputVector[which(outputVector==0)]<-0.1
  outputVector<-outputVector/nTrees
  outputVector<-log(outputVector)
  lnL<-sum(outputVector)
  return(lnL)
}