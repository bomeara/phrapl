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
Genealogies were simulated under a coalescent isolation only model where
the true history is where populations A and B coalesce after 4Ne generations, population 
A-B coalesces with C after 10Ne generations and population A-B-C coalesces with 
D after 20Ne generations. 

To access this dataset, simply load the file into R; e.g., data(TestData).  

Five objects are included:  
1. migrationArray, which lists specifications for six models to be analyzed 
(model number 1 corresponds to the correct model); 
2. migrationArrayMap, which lists parameter indexes for each model included 
in migrationArray; 
3. trees, which contains the 10 full gene trees in newick format; 
4. assignFile, which gives assignments of samples to populations;
5. popAssignments, which gives a list of vectors describing how the dataset 
should be subsampled (e.g., list(c(3,3,3)) means that one 
    subsampling procedure should be carried out, where 3 individuals are sampled
    from populations A, B, and C).

Subsampling (using PrepSubsampling) and weights calculation 
(using GetPermutationWeightsAcrossSubsamples) can first be carried out 
using assignFile popAssignments, and trees.

A phrapl GridSearch can be run using the output from these functions in 
addition to migrationArray and migrationArrayMap.
}
\examples{
  data(TestData)
  \dontrun{
  #Run PrepSubsampling on this dataset to get
  #10 subsample iterations per tree. Include 
  #everything in the working directory
  observedTrees<-PrepSubsampling(assignmentsGlobal=assignFile,
    observedTrees<-trees,popAssignments=popAssignments,
    subsamplesPerGene=10)

  #To get the weights
  subsampleWeights.df<-GetPermutationWeightsAcrossSubsamples(
    popAssignments=popAssignments,observedTrees=observedTrees)

  #To run a grid search
  result<-GridSearch(migrationArrayMap=migrationArrayMap,migrationArray=
    migrationArray,modelRange=1:6,popAssignments=popAssignments,
    observedTrees=observedTrees,subsampleWeights.df=subsampleWeights.df,
    subsamplesPerGene=10,totalPopVector=rep(20,10),nTrees=1e3,
    print.results=TRUE,return.all=TRUE)

  #To save the results to a file
  save(list=ls(), file="phraplOutput.rda")

  #To generate an output table
  totalData<-ConcatenateResults(rdaFiles="phraplOutput.rda",
  	outFile="AIC_table.txt")
    }
}
\keyword{example}
\keyword{TestData}
