\name{searchDiscreteModelSpaceOptim}
\alias{searchDiscreteModelSpaceOptim}
\title{
searchDiscreteModelSpaceOptim
}
\description{
Optimizes over model space using a genetic algorithm. This function could be used as an alternative to exhaustive search. The searching will be within an interval for a minimum or maximum for one-dimensional unconstrained function optimization.
}
\usage{
searchDiscreteModelSpaceOptim(migrationArrayMap, migrationArray, popVector, print.ms.string = FALSE, badAIC = 1e+14, maxParameterValue = 100, nTrees = 1, msLocation = "/usr/local/bin/ms", compareLocation = "comparecladespipe.pl", assign = "assign.txt", observed = "observed.txt", unresolvedTest = TRUE, debug = FALSE, method = "BFGS", itnmax = NULL, pop.size = 50, print.results = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{migrationArrayMap}{
%%     ~~Describe \code{migrationArrayMap} here~~
}
  \item{migrationArray}{
%%     ~~Describe \code{migrationArray} here~~
}
  \item{popVector}{
%%     ~~Describe \code{popVector} here~~
}
  \item{print.ms.string}{
%%     ~~Describe \code{print.ms.string} here~~
}
  \item{badAIC}{
%%     ~~Describe \code{badAIC} here~~
}
  \item{maxParameterValue}{
%%     ~~Describe \code{maxParameterValue} here~~
}
  \item{nTrees}{
%%     ~~Describe \code{nTrees} here~~
}
  \item{msLocation}{
%%     ~~Describe \code{msLocation} here~~
}
  \item{compareLocation}{
%%     ~~Describe \code{compareLocation} here~~
}
  \item{assign}{
%%     ~~Describe \code{assign} here~~
}
  \item{observed}{
%%     ~~Describe \code{observed} here~~
}
  \item{unresolvedTest}{
%%     ~~Describe \code{unresolvedTest} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{itnmax}{
%%     ~~Describe \code{itnmax} here~~
}
  \item{pop.size}{
%%     ~~Describe \code{pop.size} here~~
}
  \item{print.results}{
%%     ~~Describe \code{print.results} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
O'Meara B, Carstens B, Jackson N, Morales-Garcia A. (in prep) 
}
\author{
Brian O'Meara
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (migrationArrayMap, migrationArray, popVector, print.ms.string = FALSE, 
    badAIC = 1e+14, maxParameterValue = 100, nTrees = 1, msLocation = "/usr/local/bin/ms", 
    compareLocation = "comparecladespipe.pl", assign = "assign.txt", 
    observed = "observed.txt", unresolvedTest = TRUE, debug = FALSE, 
    method = "BFGS", itnmax = NULL, pop.size = 50, print.results = FALSE, 
    ...) 
{
    Domains <- matrix(ncol = 2, nrow = 3)
    Domains[1, ] <- range(migrationArrayMap$collapseMatrix.number)
    Domains[2, ] <- range(migrationArrayMap$n0multiplierMap.number)
    Domains[3, ] <- range(migrationArrayMap$migrationArray.number)
    results <- genoud(searchContinuousModelSpaceOptim, nvars = 3, 
        max = FALSE, starting.values = c(1, 1, 1), MemoryMatrix = TRUE, 
        boundary.enforcement = 2, data.type.int = TRUE, Domains = Domains, 
        migrationArrayMap = migrationArrayMap, migrationArray = migrationArray, 
        popVector = popVector, print.ms.string = print.ms.string, 
        badAIC = badAIC, maxParameterValue = maxParameterValue, 
        nTrees = nTrees, msLocation = msLocation, compareLocation = compareLocation, 
        assign = assign, observed = observed, unresolvedTest = unresolvedTest, 
        debug = debug, method = method, itnmax = itnmax, pop.size = pop.size, 
        print.results = print.results, ...)
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
