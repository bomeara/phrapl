
#Converts individual numbers to population labels when trees are phylo
ConvertAlleleNumbersToPopulationLetters.phylo <- function(phy,assignFrame){
	phy$tip.label <- as.character(assignFrame$popLabel[match(phy$tip.label, assignFrame$indivTotal)])
}

#Converts individual numbers to population labels when trees are raw strings
ConvertAlleleNumbersToPopulationLetters.raw<-function(phy,assignFrame){
	phy<-gsub(pattern=";", replacement="", x=gsub(pattern='\\:\\d+\\.*\\d*', replacement="", 
		x= phy, perl=TRUE)) #remove branch lengths
	split.tree<-strsplit(phy,"")[[1]]
	#split.tree<-strsplit(phy, "[^\\d]", perl=TRUE)[[1]]
	#split.tree<-split.tree[which(nchar(split.tree)>0)]
	counter=1
	#Put together double digit numbers
	for(i in 1:length(split.tree)){
		if(counter<length(split.tree)){
			if(sum((grep("\\d",split.tree[counter])),(grep("\\d",split.tree[counter+1]))) == 2){
				split.tree[counter]<-paste(split.tree[counter],split.tree[counter+1],sep="")
				split.tree<-split.tree[-(counter+1)]
				counter<-counter+2
			}else{
				counter<-counter+1
			}
		}
	}
	allele.indexes<-grep("\\d",split.tree,perl=TRUE) #indexes of taxa to transform
	#split.tree <- as.character(assignFrame$popLabel[match(split.tree,assignFrame$indivTotal)])
	for(i in allele.indexes){
		split.tree[i]<-as.character(assignFrame$popLabel[match(split.tree[i],assignFrame$indivTotal)])
	}
	return(paste(split.tree,collapse=""))
}

################################UNUSED FUNCTIONS FOR REPLACING PERL WITH R##################################


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
	outputstringMS<-paste(msLocation,sprintf("%i",msCallInfo$nsam),sprintf("%i",msCallInfo$nreps),msCallInfo$opts,
		" | grep ';' > mstrees.txt", sep=" ") 
	system(outputstringMS,intern=TRUE)
	
	#Retrieve simulated trees and population assignments
	phy<-read.table("mstrees.txt",stringsAsFactors=FALSE)
	assignFrame<-CreateAssignment.df(popAssignments[[1]])

	#Convert individual numbers to population letters and export 
	for(i in 1:length(phy[,])){
		phy[i,1]<-ConvertAlleleNumbersToPopulationLetters.raw(phy=phy[i,1],assignFrame=assignFrame)
	}

	#If going to subsample these trees, convert to phylo objects
	if(length(popAssignments) > 1){
		write.table(phy,"mstrees.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,eol=";\n")
		phy.pops<-read.tree("mstrees.txt")
	}

#	phy.phylo<-read.tree("mstrees.txt") #Also load as phylo, so can drop tips later on if necessary
#	phy.phylo$tip.label<-ConvertAlleleNumbersToPopulationLetters.phylo(phy.phylo,assignFrame)

	#Get observed trees
	observed<-read.table(observed,stringsAsFactors=FALSE)

	#Calculate popAssignment indexes for the matching vector
	treesPerLocus<-subsamplesPerGene * length(popAssignments) #number of trees per locus
	nLoci<-length(observed[,]) / treesPerLocus #number of loci
	outputVector<-array(NA,length(observed[,]))
	counter1<-1
	for(n in 1:nLoci){
		counter2<-1
		for(m in 1:(length(observed[,]) / subsamplesPerGene)){
			outputVector[counter1]<-counter2
			counter1<-counter1+1	
			if(m%%subsamplesPerGene==0){
				counter2<-counter2+1
			}
		}
	}

	#Calculate number of matches for each observed tree to each simulated tree
	currentPhy<-phy
	for(q in 1:length(popAssignments)){ #for each popAssignment size class
		#get clades for each simulated tree for this size class and convert indiv numbers to pop letters
		#Don't need to convert numbers to letter here as it was done above
		cladesAllPhy<-lapply(currentPhy[,],GetCladeNamesFromNewickString,assignFrame=assignFrame)
		for(r in which(outputVector==q)){ #for each observed subsample within a size class
			matches<-0
			cladesGene<-GetCladeNamesFromNewickString(observed[r,1],assignFrame=assignFrame,
				convertIndivNumsToPopLetters=TRUE,getOutDegreeOfPolytomies=TRUE) #get clades for current obs tree
			for(s in 1:length(currentPhy[,])){ #for each simulated tree
				matches<-matches + GetScoreOfSingleTree(cladesMS=cladesAllPhy[[s]], phyMS=currentPhy[s,1],
					cladesGene=cladesGene[[1]], phyGene=observed[r,1],polytomyDegreesForCorrection=cladesGene[[2]])
			}
			outputVector[r]<-matches
		}		
		##When we reach the next size class, drop tips from simulated trees (to do this, need to convert trees to phylo)
		if(q < length(popAssignments)){ 
			currentPhy<-phy.pops
			for(t in 1:length(currentPhy)){			
				currentPhy[[t]]<-SubsampleMSTree(phy=currentPhy[[t]],popVectorOrig=popAssignments[[q]],
					popVectorFinal=popAssignments[[q+1]])
			}
			write.tree(currentPhy,file="mstrees.txt")
			currentPhy<-read.table("mstrees.txt",stringsAsFactors=FALSE)
			assignFrame<-CreateAssignment.df(popAssignments[[q+1]])
		}
	}
	
	return(outputVector)

}

#This following set of 3 functions imports and handles trees as newick strings in attempt to speed up the 
#process of getting weights from subsamples.This attempt has resulted in failure: this is actually slower
#than using ape phylo objects to accomplish this (see the non-raw analogue functions above).
#When the armageddon comes and ape is no longer available, these functions will be there to fill the void.

#Get weights for each subsample based on the number of matches per permutation
GetPermutationWeightsAcrossSubsamples.raw<-function(phy="observed.tre",popAssignments,subsamplesPerGene,subsamplePath){	
	phy<-read.table(paste(subsamplePath,phy,sep=""),stringsAsFactors=FALSE)
	treesPerLocus<-subsamplesPerGene * length(popAssignments)
	nLoci<-nrow(phy) / treesPerLocus
	treeCounter<-1
	subsampleWeights<-c()
	
	for(y in 1:nLoci){
		subsampleSizeCounter<-1
		for(z in 1:(nrow(phy) / nLoci)){
			subsampleWeights<-append(subsampleWeights,GetPermutationWeights(phy=phy[treeCounter,],
				popVector=popAssignments[[subsampleSizeCounter]],subsamplePath=subsamplePath))
			treeCounter<-treeCounter + 1
			if(z%%subsamplesPerGene == 0){ #increase to next popAssignments size class
				subsampleSizeCounter<-subsampleSizeCounter + 1
			}
		}
	}
	return(subsampleWeights)
}

#Change labels for a population in a single newick tree to mach a new label permutation
GetPermutatedLabels.raw<-function(phy,assignFrame,popLabels=popLabels[1],newLabels=try.label){
	phy.split<-strsplit(phy,"")[[1]]
	taxonLabels<-grep("\\w",phy.split,value=TRUE)
	oldLabels<-taxonLabels[order(as.numeric(taxonLabels))][which(assignFrame[,1] == popLabels[1])]
	oldLabelIndex<-array(NA,dim=length(oldLabels))
	for(i in 1:length(oldLabels)){
		oldLabelIndex[i]<-which(phy.split==oldLabels[i])
	}
	for(i in 1:length(oldLabels)){
		phy.split[oldLabelIndex[i]]<-newLabels[i]
	}
	return(paste(phy.split,collapse=""))
}

#Get weights for a given tree based on the number of matches per permutation
GetPermutationWeights.raw<-function(phy,popVector,subsamplePath){
	phy<-gsub(pattern=";",replacement="",x=gsub(pattern='\\:\\d+\\.*\\d*',replacement="", 
		x=phy,perl=TRUE)) #remove branch lengths
	assignFrame<-CreateAssignment.df(popVector)
	popLabels<-unique(assignFrame[,1])
	#Make list of label permutations for each population
	tips.permute<-list()
	for(i in 1:length(popLabels)){
		tips.permute[[length(tips.permute) + 1]]<-as.matrix(perms(length(assignFrame[,2][which(assignFrame[,1] == 
			popLabels[i])])))
		tips.permute[[i]]<-tips.permute[[i]] + (((assignFrame[,2][which(assignFrame[,1] == popLabels[i])])[1]) - 1)
	}
	numberOfPermutations<-0
	phyOriginal<-phy
	cladesMS<-GetCladeNamesFromNewickString(phy,assignFrame)
	matches<-0
	if(length(popLabels) > 10){
		stop("This function can only handle 10 populations")
	}
	
	#Start permutation of populations	
	for(j in 1:ncol(tips.permute[[1]])){
		phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
			popLabels=popLabels[1],newLabels=tips.permute[[1]][,j])
					
		if(length(popLabels) > 1){	
			for(k in 1:ncol(tips.permute[[2]])){
				phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
					popLabels=popLabels[2],newLabels=tips.permute[[2]][,k])			

				if(length(popLabels) > 2){	
					for(l in 1:ncol(tips.permute[[3]])){
						phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
							popLabels=popLabels[3],newLabels=tips.permute[[3]][,l])	

						if(length(popLabels) > 3){
							for(m in 1:ncol(tips.permute[[4]])){
								phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
									popLabels=popLabels[4],newLabels=tips.permute[[4]][,m])
							
								if(length(popLabels) > 4){	
									for(n in 1:ncol(tips.permute[[5]])){
										phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
											popLabels=popLabels[5],newLabels=tips.permute[[5]][,n])				
									
										if(length(popLabels) > 5){	
											for(o in 1:ncol(tips.permute[[6]])){
												phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
													popLabels=popLabels[6],newLabels=tips.permute[[6]][,o])	
											
												if(length(popLabels) > 6){	
													for(p in 1:ncol(tips.permute[[7]])){
														phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
															popLabels=popLabels[7],newLabels=tips.permute[[7]][,p])				
														
														if(length(popLabels) > 7){	
															for(q in 1:ncol(tips.permute[[8]])){
																phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
																	popLabels=popLabels[8],newLabels=tips.permute[[8]][,q])	
																
																if(length(popLabels) > 8){
																	for(r in 1:ncol(tips.permute[[9]])){
																		phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
																			popLabels=popLabels[9],newLabels=tips.permute[[9]][,r])
																		
																		if(length(popLabels) > 9){	
																			for(s in 1:ncol(tips.permute[[10]])){
																				phy<-GetPermutatedLabels(phy=phy,assignFrame=assignFrame,
																					popLabels=popLabels[10],newLabels=tips.permute[[10]][,s])				

																				if(length(popLabels) == 10){
																					cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
																					matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
																					numberOfPermutations<-numberOfPermutations + 1
																				}
																			}
																		}
																		if(length(popLabels) == 9){
																			cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
																			matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
																			numberOfPermutations<-numberOfPermutations + 1
																		}
																	}
																}
																if(length(popLabels) == 8){
																	cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
																	matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
																	numberOfPermutations<-numberOfPermutations + 1
																}
															}
														}
														if(length(popLabels) == 7){
															cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
															matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
															numberOfPermutations<-numberOfPermutations + 1
														}
													}
												}
												if(length(popLabels) == 6){
													cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
													matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
													numberOfPermutations<-numberOfPermutations + 1
												}
											}
										}
										if(length(popLabels) == 5){
											cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
											matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
											numberOfPermutations<-numberOfPermutations + 1
										}
									}
								}
								if(length(popLabels) == 4){
									cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
									matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
									numberOfPermutations<-numberOfPermutations + 1
								}
							}
						}
						if(length(popLabels) == 3){
							cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
							matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
							numberOfPermutations<-numberOfPermutations + 1
						}
					}
				}
				if(length(popLabels) == 2){
					cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
					matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
					numberOfPermutations<-numberOfPermutations + 1
				}
			}
		}
		if(length(popLabels) == 1){
			cladesGene<-GetCladeNamesFromNewickString(phy,assignFrame)
			matches<-matches + GetScoreOfSingleTree(cladesMS,phyOriginal,cladesGene,phy)
			numberOfPermutations<-numberOfPermutations + 1
		}
	}
	weight<-matches / numberOfPermutations
	write.table(weight,file=paste(subsamplePath,"subsampleWeights.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,
				col.names=FALSE,append=TRUE)
	return(weight)
}	


#Add specific descriptor columns to results (specific for simulation checking)
AddDescriptors<-function(dataframe=totalData){           
	#Get treatment descriptor columns
   	whichdataframe<-as.matrix(rep(datasetsVec[dataset],length(nrow(dataframe)))) #make dataframe descriptor column
   	migTreat<-as.matrix(rep(migVec[dataset],length(nrow(dataframe)))) #make migration column
   	divTreat<-as.matrix(rep(divVec[divRep],length(nrow(dataframe)))) #make divergence column
   	treeTreat<-as.matrix(rep(treeVec[treeRep],length(nrow(dataframe)))) #make tree number column
   	subnumTreat<-as.matrix(rep(subsampleNumVec[subsampleNumRep],length(nrow(dataframe)))) #make simulation rep column
   	subsampleTreat<-as.matrix(rep(subsampleRep,length(nrow(dataframe)))) #make subsample rep column
	
	#Combine dataframe across models with treatment descriptor columns
	dataframe<-cbind(whichdataframe,migTreat,divTreat,treeTreat,subnumTreat,subsampleTreat,dataframe)
	colnames(dataframe)<-c("dataset","migration","divTime","nTrees","subNum","subsample",colnames(dataframe[,7:ncol(dataframe)]))
    return(dataframe)
}   