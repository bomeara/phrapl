#ported from perl
library(ape)

#left to do
#permutations? How to correct for (A,(A,(A,A))) matching (1,(2,(3,4))) and (1,(2,(4,3)))
#decide on how to do subsamples for MS in order to extrapolate to final AIC
#add extrapolation code to grid search
#test


ConvertAlleleNumbersToPopulationLetters <- function(phy, popVector) {
	assignFrame <- CreateAssignment.df(popVector)
	phy$tip.label <- assignFrame$popLabel[match(phy$tip.label, assignFrame$indivTotal)]
}

#assumes same number of populations in the two popVectors, just that we subsample from one to the other
#though it would be possible to have popVectorFinal have zero entries, thus effectively deleting that population
SubsampleMSTree <- function(phy, popVectorOrig, popVectorFinal) {
	taxaToDrop <- c()
	minSample <- 1
	maxSample <- 0
	for (population in sequence(length(popVectorOrig))) {
		maxSample <- minSample + popVectorOrig[population] - 1
		individualsThisPop <- c(minSample:maxSample)
		sampleSizeDiff <- popVectorOrig[population] - popVectorFinal[population]
		if (sampleSizeDiff > 0) {
			taxaToDrop <- append(taxaToDrop, individualsThisPop[-c(1: sampleSizeDiff)])
		}
		minSample <- maxSample + 1
	}
	phy <- drop.tip(phy, taxaToDrop)
	phy$tip.label <- order(order(as.numeric(phy$tip.label))) #we have a set of taxon labels like 3, 4, 6, 7, 9, 10 and convert these to 1, 2, 3, 4, 5, 6: keeping the order (since this corresponds to the assignment to population) but filling in gaps
	return(phy)
}

#gets a list of all clades in the tree. For the taxa descended from each clade, sorts alphabetically and then makes them a string
GetAndBindLabel <- function(phy) { 
	#note the sorting here
	return( paste( sort( x$tip.label ), collapse="_" ) )
}

GetOutDegreeOfPolytomies <- function(phy) {
	descendantCounts <- table(phy$edge[,1])
	descendantCounts <- unname(descendantCounts[whichDescendantCounts>2])
	return(descendantCounts)
}

#assumes:
#1. You have already run ConvertAlleleNumbersToPopulationLetters so these trees have letters
#2. You have already made them have the same size (do SubsampleMSTree if needed)
GetScoreOfSingleTree <- function(phyMS, phyGene) {
	cladesMS <- simplify2array(sapply(subtrees(phyMS), GetAndBindLabel))
	cladesGene <- simplify2array(sapply(subtrees(phyGene), GetAndBindLabel))
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


#DO NOT USE THIS YET; IT HAS YET TO BE TESTED, AND DOESN'T DEAL WITH
#ASSIGNMENTS YET


does.agree.unrooted<-function(phyA,phyB) { #so you need to have added outgroup if you want this rooted
	#break down trees into their splits. If they are compatible, the max of the compatibility matrix will be 0.
	#they must have the same tips
	if(max(compatible(as.splits(prop.part(phyA,phyB))))==1) {
		return(FALSE)
	}
	else {
		return(TRUE)
	}
}

unresolvedScaleFactor<-function(phy) {
	numberUnresolvedEdges<-(Ntip(phy)-1) - Nnode(phy)
	return(1.0/(3^numberUnresolvedEdges))
}

getScoreOfPair<-function(simulatedTree, observedTree) {
	if (does.agree.unrooted(observedTree,simulatedTree)) {
		return(unresolvedScaleFactor(observedTree))
	}
	else {
		return(0)
	}
}

getScoresOfManySimulated<-function(simulatedTreesVector,observedTree) {
  if (class(simulatedTreesVector)=="multiPhylo") {
     simulatedTreesVector<-unclass(simulatedTreesVector)
  }
	return(sapply(simulatedTreesVector,getScoreOfPair, observedTree=observedTree ,simplify=TRUE))
}

addOutgroup<-function(phy,outgroupName="outgroup") {
  return(bind.tree(read.tree(text=paste('(',outgroupName,',',outgroupName,');',sep="")),phy,1))
}

pruneToSameLeafSet<-function(simulatedTree,observedTree) {
  leavesSim<-simulatedTree$tip.label
  leavesToPrune<-leavesSim[which(is.na( match( leavesSim, observedTree$tip.label)))]
  if(length(leavesToPrune)>0) {
    print(leavesToPrune)
    simulatedTree<-drop.tip(simulatedTree,tip=leavesToPrune)
  }
  return(simulatedTree)
}

pruneAndOutgroup<-function(simulatedTree,observedTree,outgroupName="outgroup") {
  simulatedTree<-pruneToSameLeafSet(simulatedTree,observedTree)
  simulatedTree<-addOutgroup(simulatedTree,outgroupName=outgroupName)
  return(simulatedTree)
}

pruneAndOutgroupForManySimulated<-function(simulatedTreesVector,observedTree) {
  if (class(simulatedTreesVector)=="multiPhylo") {
     simulatedTreesVector<-unclass(simulatedTreesVector)
  }
  return(lapply(simulatedTreesVector,pruneAndOutgroup, observedTree=observedTree))
}
