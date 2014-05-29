\name{GridSearch}
\alias{GridSearch}
\alias{Gridsearch}
\alias{gridsearch}
\alias{grid.search}
\title{
Computes AIC for a given model or models on a fixed set of parameters
}
\description{
This function takes a given model or models and a grid of parameters and returns
the AIC at each parameter. In testing, it was found to work better than 
traditional optimization approaches. However, it can still do a heuristic 
search.
}
\usage{
GridSearch(modelRange = c(1:length(migrationArray)), migrationArrayMap, 
migrationArray, popAssignments, badAIC = 1e+14, 
maxParameterValue = 100, nTrees = 1e5, msPath = "ms", 
comparePath = system.file("extdata","comparecladespipe.pl", package = "phrapl"),
observedPath = "observed.tre", subsampleWeightsPath = "subsampleWeights.txt", 
unresolvedTest = TRUE, print.ms.string = FALSE, print.results = FALSE, 
print.matches = FALSE, debug = FALSE, method = "nlminb", itnmax = NULL, 
ncores = 1, results.file = NULL, maxtime = 0, maxeval = 0, return.all = TRUE, 
numReps = 0, startGrid = NULL, 
collapseStarts = c(0.3, 0.58, 1.11, 2.12, 4.07, 7.81, 15), 
n0Starts = c(0.1, 0.5, 1, 2, 4), 
migrationStarts = c(0.1, 0.22, 0.46, 1, 2.15, 4.64, 10), 
subsamplesPerGene = 1, totalPopVector = totalPopVector, summaryFn = "mean", 
saveNoExtrap = FALSE, doSNPs = FALSE, ...)
}
\arguments{
   \item{modelRange}{
Integer vector: which models to examine. Do not specify to use default of all
models.
}
 \item{migrationArrayMap}{
A data.frame containing information about all the models
}
  \item{migrationArray}{
List containing all the models
}
  \item{popAssignments}{
A list of vectors (typically only one vector will be specified) that define the 
number of individuals per population included in the observed tree file 
(usually these will be subsampled trees). 
Defining popAssignments as list(c(4,4,4)) for example means that there 12 tips 
per observed tree, with 4 tips per population.
}
  \item{badAIC}{
In case of failure (such as trying a parameter outside a bound), this allows 
returning of suboptimal but still finite number. Mostly used for heuristic 
searches.
}
  \item{maxParameterValue}{
A bound for the maximum value for any parameter.
}
  \item{nTrees}{
Integer: the number of trees to simulate in ms.
}
  \item{msPath}{
Path to the local installation of ms; typing this string on the command line 
should result in ms running.
}
  \item{comparePath}{
Path to the local placement of the compareCladesPipe.pl perl script, including 
that script name.
}
  \item{observedPath}{
Path to the local instance of the observed trees, including the name of the 
file.
 }
  \item{subsampleWeightsPath}{
Path to the local instance of the weights for the trees, including the name of 
the file.
}
  \item{unresolvedTest}{
Boolean: deal with unresolved gene trees by looking for partial matches and 
correcting for that.
}
  \item{print.ms.string}{
Mostly for debugging, Boolean on whether to verbosely print out the calls to ms.
}
  \item{print.results}{
Mostly for debugging, Boolean on whether to verbosely print out the results.
}
  \item{print.matches}{
Mostly for debugging, Boolean on whether to verbosely print out the matches.}
  \item{debug}{
Whether to print out additional debugging information.
}
  \item{method}{
For heuristic searches, which method to use. ?optim for more information.
}
  \item{itnmax}{
For heuristic searches, how many steps.
}
  \item{ncores}{
Allows running on multiple cores. Not implemented yet.
}
  \item{results.file}{
File name for storage of results.
}
  \item{maxtime}{
Maximum run time for heuristic search.
}
  \item{maxeval}{
Maximum number of function evaluations to run for heuristic search.
}
  \item{return.all}{
Boolean: return just the AIC scores or additional information.
}
  \item{numReps}{
For heuristic searches, number of starting points to try.
}
  \item{startGrid}{
Starting grid of parameters to try. Leave NULL to let program create this.
}
  \item{collapseStarts}{
Vector of starting values for collapse parameters.
}
  \item{n0Starts}{
Vector of starting values for n0.
}
  \item{migrationStarts}{
Vector of starting values for migration rates.
}
  \item{subsamplesPerGene}{
How many subsamples to take per gene
}
  \item{totalPopVector}{
Overall number of samples in each population before subsampling.
}
  \item{summaryFn}{
Way to summarize results across subsamples.
}
  \item{saveNoExtrap}{
Boolean to tell whether to save extrapolated values. FALSE by default.
}
  \item{doSNPs}{
Boolean to tell whether to use the SNPs model: count a single matching edge on
the simulated tree as a full match. FALSE by default.
}
  \item{\dots}{
Other items to pass to heuristic search functions.
}
}
\details{
We recommend using the grid, not the heuristic search. If you are using the
SNPs model, you should have uniform weights for the gene trees unless
there is some reason you want to weight them differently.
}
\value{
If return.all==FALSE, just a vector of AIC values. Otherwise, a list with 
parameters used for the grid and AIC for each, if using grid search.
}
\author{
Brian O'Meara & Nathan Jackson
}
\note{
For more information, please see the user manual.
}
\keyword{ ~grid }
\keyword{ ~gridsearch }
