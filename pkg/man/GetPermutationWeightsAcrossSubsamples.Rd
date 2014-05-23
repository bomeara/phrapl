\name{GetPermutationWeightsAcrossSubsamples}
\alias{GetPermutationWeightsAcrossSubsamples}
\title{
Calculate tip label degeneracy weights for a set of trees
}
\description{
This function calculates weights for a tree or set of trees (typically subsampled trees) based on
the degeneracy of tip lables (i.e., the proportion of times that permuting labels across a
tree results in the same topology). These weights are applied during a grid or heuristic search in phrapl 
to obtain a more accurate estimate of the log likelihood of an observed tree given trees expected under a given model.
}
\usage{
GetPermutationWeightsAcrossSubsamples(popAssignments,subsamplePath="./",inputFile="observed.tre",
outputFile="subsampleWeights.txt")
}
\arguments{
  \item{popAssignments}{
A list of vectors (typically only one vector will be specified) that define the number of tips per 
population included in the observed tree file (usually these will be subsampled trees).
}
  \item{subsamplePath}{
Path to the directory from which input files are accessed and to which output files are printed.
}
  \item{inputFile}{
Name of the file that includes the subsampled trees
}
  \item{outputFile}{
Name of the file to which the weights will be printed (its length will equal the number of subsampled trees)
}
\note{
Because this function exhaustively permutes all possible tip lables across a tree, it can take a long time
for large subsampled trees. For more information, please see the user manual.
}
\author{
Nathan Jackson & Brian O'Meara
}
\note{
For more information, please see the user manual.
}
\keyword{ ~weights }
\keyword{ ~permutation }