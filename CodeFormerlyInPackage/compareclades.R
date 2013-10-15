#ported from perl
library(ape)
library(phangorn)

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
