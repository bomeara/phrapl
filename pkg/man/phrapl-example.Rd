\name{phrapl-example}
\alias{phrapl-example}
\alias{example}
\alias{TestData}
\alias{migrationArray}
\alias{migrationArrayMap}
\alias{trees}
\alias{assignFile}
\alias{popAssignments}
\alias{subsampledTrees}
\alias{subsampleWeights}
\title{
A description for applying phrapl to a test dataset 
}
\description{
A test dataset called TestData.rda can be found in the data directory included 
in the phrapl package. This dataset comprises 
10 gene trees containing three populations ("A","B", and "C"), with 20 samples 
per population, and a single outgroup sample ("D"). 
Genealogies were simulated under a coalescent isolation-with-migration model in
which populations A and B share alleles at a rate of 
Nm = 0.5 and coalesce after Ne generations. Population A-B coalesces with C
after 2.5Ne generations and population A-B-C coalesces with 
D after 5Ne generations (no migration involved). 

To access this dataset, simply load the file into R; e.g., data(TestData).  

Seven objects are included:  
1. migrationArray, which lists specifications for six models to be analyzed 
(model number 2 corresponds to the correct model); 
2. migrationArrayMap, which lists parameter indexes for each model included 
in migrationArray; 
3. trees, which contains the 10 full gene trees in newick format; 
4. assignFile, which gives assignments of samples to populations;
5. popAssignments, which gives a list of vectors describing how the dataset 
should be subsampled (e.g., list(c(4,4,4)) means that one 
    subsampling procedure should be carried out, where 4 individuals are sampled
    from populations A, B, and C);
6. subsampledTrees, which contains subsampled gene trees (10 subsampled trees 
for each original gene tree, sampled at 4 samples per population;
	the outgroup has been pruned in each tree);
7. subsampleWeights, which gives a weight for each subsample topology that is 
used in phrapl to calculate log likelihoods.

A phrapl GridSearch can be run using migrationArray, migrationArrayMap, 
popAssignments, subsampledTrees, and subsampleWeights. 
To perform subsampling (PrepSubsampling) and weights calculation 
(GetPermutationWeightsAcrossSubsamples) ahead of doing the grid search, 
use assignFile and trees.

Phrapl functions input assignFile, trees, subsampledTrees, and subsampleWeights 
as files rather than as R objects. Thus, to use these four 
objects, first export them to a file (note you will need to install and load the
ape package to export the trees).

}
\examples{
  data(TestData)
  \dontrun{
  #load ape
  library(ape)
  
  #Write trees to file
  write.tree(trees,file="trees.tre")

  #Write subsampled trees to file
  write.tree(subsampledTrees,file="observed.tre")

  #Write assignFile to file
  write.table(assignFile,file="cladeAssignments.txt",
    quote=FALSE,row.names=FALSE)

  #Write subsampleWeights to file
  write.table(subsampleWeights,file=
    "subsampleWeights.txt",quote=FALSE,row.names=FALSE,
    col.names=FALSE)

  #To run PrepSubsampling on this dataset, include 
  #everything in the working directory
  PrepSubsampling(popAssignments=popAssignments,
    subsamplesPerGene=10)

  #To get the weights
  GetPermutationWeightsAcrossSubsamples(popAssignments=
    popAssignments)

  #To run a grid search
  result<-GridSearch(migrationArrayMap=
    migrationArrayMap,migrationArray=
    migrationArray,modelRange=1:6,
    popAssignments=popAssignments,subsamplesPerGene=
    10,nTrees=1e4,print.results=TRUE,return.all=TRUE)

  #To save the results to a file
  save(list=ls(), file="phraplOutput.rda")

  #To generate an output table
  totalData<-ConcatenateResults(outFile=
    "AIC_table.txt")
    }
}
\keyword{example}
\keyword{TestData}
