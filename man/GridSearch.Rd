\name{GridSearch}
\alias{GridSearch}
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
GridSearch(modelRange = c(1:length(migrationArray)), migrationArrayMap=NULL,
migrationArray, popAssignments, badAIC = 1e+14,
maxParameterValue = 20, nTrees = 1e5, msPath = system.file("msdir", "ms", package = "P2C2M"),
comparePath = system.file("extdata","comparecladespipe.pl", package = "phrapl"),
observedTrees, subsampleWeights.df = NULL, doWeights = TRUE,
unresolvedTest = TRUE, print.ms.string = FALSE, print.results = TRUE,
print.matches = FALSE, debug = FALSE, method = "nlminb", itnmax = NULL,
ncores = 1, results.file = NULL, maxtime = 0, maxeval = 0, return.all = TRUE,
numReps = 0, startGrid = NULL,
collapseStarts = c(0.3, 0.58, 1.11, 2.12, 4.07, 7.81, 15),
n0Starts = c(0.1, 0.5, 1, 2, 4),
growthStarts=c(0.30,0.58,1.11,2.12,4.07,7.81,15.00),
migrationStarts = c(0.1, 0.22, 0.46, 1, 2.15, 4.64),
subsamplesPerGene = 1, totalPopVector = NULL, summaryFn = "mean",
saveNoExtrap = FALSE, doSNPs = FALSE, nEq=100, setCollapseZero=NULL,
dAIC.cutoff=2, rm.n0 = TRUE, popScaling = NULL, checkpointFile = NULL,
parameter.ambiguous = FALSE, addedEventTime = NULL, addedEventTimeAsScalar = TRUE,
optimization = "grid", genoudPopSize = 25, numGridStartVals = 25,
solutionTolerance = 1, skipGrid = FALSE, usePhyclust=FALSE, ...)
}
\arguments{
   \item{modelRange}{
Integer vector: which models to examine. Do not specify to use default of all
models.
}
 \item{migrationArrayMap}{
A data.frame containing information about all the models. Only required for
heuristic search, not for grid search.
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
  \item{observedTrees}{
Multiphylo object of the empirical trees.
 }
  \item{subsampleWeights.df}{
A dataframe of the weights for each subsample. If this is NULL, it is computed
within GridSearch.
}
  \item{doWeights}{
In no subsampleWeights.df object is called and doWeights = TRUE, subsample weights
will be calculated for each tree prior to AIC calculation.
}
  \item{unresolvedTest}{
Boolean: deal with unresolved gene trees by looking for partial matches and
correcting for that.
}
  \item{print.ms.string}{
Mostly for debugging, Boolean on whether to verbosely print out the calls to ms.
}
  \item{print.results}{
If TRUE, after each simulation cycle, parameters and AICs are printed to the screen.
}
  \item{print.matches}{
Mostly for debugging, Boolean on whether to verbosely print out the matches.
}
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
  \item{growthStarts}{
Vector of starting values for growth parameters.
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
  \item{nEq}{
If no simulated trees match the observed trees, the frequentist estimate of
the matching proportion is exactly zero (flip a coin 5 times, see no heads, so
estimate probability of heads is zero). This would have an extreme effect: if
any gene trees don't match, the model has no likelihood. A better approach is
to realize that a finite set of samples gives finite information. We assume
that the probability of a match, absent data, is 1/number of possible gene
trees, and this is combined with the empirical estimate to give an estimate of
the likelihood. A question is how much weight to put on this pre-existing
estimate, and that is set by nEq. With its default of 100, very low weight is
placed on this: it's equivalent in the info present in nTrees=100, and since
the actual nTrees is 10000 or more, it has very little impact.
}
  \item{setCollapseZero}{
A vector of collapse parameters that will be set to zero (e.g., c(1,2) will set
both the first and second collapse parameters to zero). K will be adjusted
automatically to account for the specified fixed parameters.
}
  \item{dAIC.cutoff}{
A value specifying how optimized parameter values should be selected. Parameter
estimates are calculated by taking the mean parameter value across all values
within an AIC distance of dAIC.cutoff relative to the lowest AIC value. The
default is 2 AIC points. Note that this argument only applies when calculating
parameters with the ambiguous name format that was used prior to
December 2015.
}
  \item{rm.n0}{
A boolean indicating whether n0multiplier parameter estimates should be outputted
with the other estimates, even if the value is always the same.
This de-clutters output if one is not analyzing models that vary
population sizes. Note that phrapl checks whether there are multiple n0multiplier values
in the inputted migrationArray; if there are, rm.n0 is automatically switched to be FALSE.
}
  \item{popScaling}{
A vector whose length is equal to the number of loci analyzed that gives the relative
scaling of effective population size for each locus (e.g., diploid nuclear locus = 1,
X-linked locus = 0.75, mtDNA locus = 0.25). Default is equal scaling for all loci.
}
  \item{checkpointFile}{
Results can be printed to a file which acts as a checkpoint. If a job is stopped during
a GridSearch, rather than starting the analysis over, if a file was previously specified
using this argument, the search will resume at the most recent iteration printed to the
file. Currently, this argument can only be used when running a single model in GridSearch.
However, if GridSearch is taking a long time when running multiple models at once, one
can save results more often by reducing the number of models run at a time. Thus, the
checkpointFile argument is only really necessary when you've got a single model running
that can't be broken up any further.
}
  \item{parameter.ambiguous}{
If true, calculates parameters with the ambiguous name format that was used prior to
December 2015, rather than using the unambiguous parameter names format.
}
  \item{addedEventTime}{
This takes either a single value or vector of values giving absolute or relative times at which
non-coalescence events that have been added to a migrationArray occur (such events are added using
the function \code{\link{AddEventToMigrationArray}}). There should be as many values listed as
there are added events in the collapseMatrix (sorted from most recent to most ancient).
Time values can be specifed as either an absolute time period (in units of 4Ne) or as a relative time
period. When the latter is desired, the \code{addedEventTimeAsScalar} argument (see below)
must be set to \code{TRUE}, which will cause the values specified by \code{addedEventTime} to be
treated as scalars, rather than as absolute. Scalar values must be greater than zero
and less than one, and will be multiplied by whatever the next collapse time estimate happens to be
(which will typically vary across models and across a search grid). Note that the specified timing of
these events must fit the modeled chronology of coalescent and non-coalescence events in the
collapseMatrix. That is, the timing for an added non-coalescence event given by \code{addedEventTime}
must be later than the collapse event directly preceding the added event
and must be earlier than the collapse event directly following it. For this reason, it is recommended
that \code{addedEventTimes} be entered as scalars, as this guarantees that the time will always be
earlier than the collapse event that follows (but not that the time will occur after the previous collapse
event). Absolute values can also be used, but remember that in a grid search, different values
are explored, so many parameter combinations in the grid may not fit the specified event times. However,
PHRAPL filters out all parameter combinations in the parameter grid that do not adhere to the specified
chronology of events in the model. Thus, setting absolute dates may reduce the parameter space explored
for a model, but this parameter space will always fit the model.
}
  \item{addedEventTimeAsScalar}{
Dictates whether added event times set by addedEventTime are treated as scalars (when \code{TRUE}) or as
absolute values (when \code{FALSE}).
}
  \item{optimization}{
Dictates how parameters are estimated. Options are \code{"grid"}, \code{"nloptr"}, or \code{"genoud"}.
If \code{"nloptr"} or \code{"genoud"} are selected, this causes PHRAPL to first run a grid search to
obtain starting values; then further optimization occurs using an algorithm.
}
  \item{genoudPopSize}{
Specifies the population size (pop.size) value used by genoud.
}
  \item{numGridStartVals}{
Specifies the number of grid value combinations used as starting values for genoud (where grid values are
sorted by AIC such that the best models are preferentially chosen). If one carries out an initial grid search to
obtain starting values for genoud (\code{skipGrid = FALSE}), then smaller \code{numGridStartVals} values will
give more weight to the best values inferred from a grid search. Using a large number for \code{numGridStartVals}
will dampen the impact of grid search optimization. Note that if \code{numGridStartVals} is less than
\code{genoudPopSize}, genoud will automatically add in addition "individuals". If \code{numGridStartVals} exceeds
the number of grid values for a given model, genoud will simply add more individuals such that \code{genoudPopSize}
is satisfied.
}
  \item{solutionTolerance}{
This gives the \code{solution.tolerance} level used by genoud
}
  \item{skipGrid}{
If set to \code{TRUE}, the initial grid search is skipped when doing genoud optimization. In this case,
starting values are determined by randomly selecting values from the grid (the number of which is determined
by \code{numGridStartVals}).
}
\item{usePhyclust}{
If TRUE, use phyclust's version of ms rather than P2C2M (should be equivalent)
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
