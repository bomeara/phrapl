\name{PrepSubsampling}
\alias{PrepSubsampling}
\alias{Prepsubsampling}
\alias{prepsubsampling}
\alias{prep.subsampling}
\title{
Subsample tips from a tree
}
\description{
This function subsamples tips from a tree according to a specified sampling 
strategy. Subsampling is ordinarily done by specifying a fixed number of 
individuals to take per population (using popAssignments). An outgroup 
population can be included to root each subsample, but then pruned from the tree
prior to analysis using phrapl.
}
\usage{
PrepSubsampling(subsamplePath="./",assignFile="cladeAssignments.txt",
treesFile="trees.tre", outputFile="observed.tre",popAssignments,
subsamplesPerGene,nIndividualsDesired=NULL,minPerPop=1,outgroup=TRUE,
outgroupPrune=TRUE)
}
\arguments{
  \item{subsamplePath}{
Path to the directory from which input files are accessed and to which output 
files are printed.
}
  \item{assignFile}{
File that includes population assignments. This should consist of two columns: 
the first lists all individuals in the tree; the second lists populations to 
which individuals belong (e.g., "A", "B", "C"...).
The file must include a header row.
}
  \item{treesFile}{
A file that includes one or more trees to be subsampled, in newick format.
}
  \item{outputFile}{
File to which subsampled trees will be printed. If there is more than one tree 
in the treesFile, all subsamples for the first tree will be printed, then all 
the subsamples for the second tree, etc.
}
  \item{popAssignments}{
A list of vectors that gives a fixed number of subsamples to take per 
population. Defining popAssignments as list(c(4,4,4)) for example will result in
4 subsamples being taken from each of three populations. The number of 
indivduals to subsample can vary among populations. If more than one subsampling
vector is included,all subsampling iterations will first be carried out 
according to the first vector, then subsequent vectors, for a given input tree. 
}
  \item{subsamplesPerGene}{
Gives the number of subsampling iterations per tree.
}
  \item{nIndividualsDesired}{
The total number of individuals to be included in a subsample. This must only be
specified if one wishes to randomly subsample tips by only fixing the total
number to be subsampled, but not fixing the specific number to be taken per 
population (which is done by specifying popAssignments). This option must be 
specified in conjuction with minPerPop. 
}
  \item{minPerPop}{
If one wishes to subsample from a tree, but not fix the number of samples to be 
taken per population, this gives the minimum number of individuals to include 
from any population in a subsample. This is only used in conjunction with 
nIndividualsDesired. If popAssignments is defined, minPerPop will be ignored.
}
  \item{outgroup}{
Whether an outgroup is included in the tree to be used to root the subsampled 
trees. If TRUE, one sample from the outgroup population will be included in each
subsampled tree. The outgroup population is identified as the last population
included in the assignFile. 
}
  \item{outgroupPrune}{
Whether a subsampled outgroup should be pruned from the subsampled tree prior to
being printed to the output file. This is to be used in conjunction with the 
outgroup option.
 }
 }
\author{
Nathan Jackson & Brian O'Meara
}
\note{
For more information, please see the user manual.
}
\keyword{ ~subsampling }
