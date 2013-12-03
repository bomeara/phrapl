\name{exhaustiveSearchOptim}
\alias{exhaustiveSearchOptim}
\title{
exhaustiveSearchOptim
}
\description{
This function will perform an exhaustive search at every model across all model space within an interval for a minimum or maximum for one-dimensional unconstrained function optimization.
}
\usage{
exhaustiveSearchOptim(migrationArrayMap, migrationArray, popVector, badAIC = 1e+14, maxParameterValue = 100, nTrees = 1, msLocation = "/usr/local/bin/ms", compareLocation = "comparecladespipe.pl", assign = "assign.txt", observed = "observed.txt", unresolvedTest = TRUE, print.ms.string = FALSE, print.results = FALSE, debug = FALSE, method = "nlminb", itnmax = NULL, ncores = ncores, results.file = NULL, return.all = TRUE, ...)
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
  \item{print.ms.string}{
%%     ~~Describe \code{print.ms.string} here~~
}
  \item{print.results}{
%%     ~~Describe \code{print.results} here~~
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
  \item{ncores}{
%%     ~~Describe \code{ncores} here~~
}
  \item{results.file}{
%%     ~~Describe \code{results.file} here~~
}
  \item{return.all}{
%%     ~~Describe \code{return.all} here~~
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
function (migrationArrayMap, migrationArray, popVector, badAIC = 1e+14, 
    maxParameterValue = 100, nTrees = 1, msLocation = "/usr/local/bin/ms", 
    compareLocation = "comparecladespipe.pl", assign = "assign.txt", 
    observed = "observed.txt", unresolvedTest = TRUE, print.ms.string = FALSE, 
    print.results = FALSE, debug = FALSE, method = "nlminb", 
    itnmax = NULL, ncores = ncores, results.file = NULL, return.all = TRUE, 
    ...) 
{
    AIC.values <- rep(NA, length(migrationArray))
    for (i in sequence(length(migrationArray))) {
        p <- c(migrationArrayMap$collapseMatrix.number[i], migrationArrayMap$n0multiplierMap.number[i], 
            migrationArrayMap$migrationArray.number[i])
        try(AIC.values[i] <- searchContinuousModelSpaceOptim(p, 
            migrationArrayMap, migrationArray, popVector, badAIC = badAIC, 
            maxParameterValue = maxParameterValue, nTrees = nTrees, 
            msLocation = msLocation, compareLocation = compareLocation, 
            assign = assign, observed = observed, unresolvedTest = unresolvedTest, 
            print.ms.string = print.ms.string, print.results = print.results, 
            debug = debug, method = method, itnmax = itnmax, 
            ...))
        print(c(i, length(migrationArray), i/length(migrationArray), 
            AIC.values[i]))
        if (!is.null(results.file)) {
            save(list = ls(), file = results.file)
        }
    }
    return(AIC.values)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
