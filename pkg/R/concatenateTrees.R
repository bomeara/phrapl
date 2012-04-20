setwd("/Users/bomeara/Desktop/treesreformated")
library(ape)
treefiles<-system("ls -1 *.tre",intern=TRUE)
observedTrees<-list
for (i in sequence(length(treefiles))) {
  print(treefiles[i])
  #observedTrees<-append(observedTrees,read.tree(file=treefiles[i]))
  newick<-read.tree(file=treefiles[i])
  nexus<-read.nexus(file=treefiles[i])
  print(newick)
  print(nexus)
  
}