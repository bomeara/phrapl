\name{searchContinuousModelSpaceNLoptr}
\alias{searchContinuousModelSpaceNLoptr}
\title{
searchContinuousModelSpaceNLoptr
}
\description{
Optimizes the continuous parameters for a given model using nonlinear optimization algorithms.
}
\usage{
searchContinuousModelSpaceNLoptr(p, migrationArrayMap, migrationArray, popVector, badAIC = 1e+14, maxParameterValue = 100, nTrees = 1, msLocation = "/usr/local/bin/ms", compareLocation = "comparecladespipe.pl", assign = "assign.txt", observed = "observed.txt", unresolvedTest = TRUE, print.ms.string = FALSE, print.results = FALSE, debug = FALSE, method = "nlminb", itnmax = NULL, return.all = FALSE, maxtime = 0, maxeval = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
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
  \item{return.all}{
%%     ~~Describe \code{return.all} here~~
}
  \item{maxtime}{
%%     ~~Describe \code{maxtime} here~~
}
  \item{maxeval}{
%%     ~~Describe \code{maxeval} here~~
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
function (p, migrationArrayMap, migrationArray, popVector, badAIC = 1e+14, 
    maxParameterValue = 100, nTrees = 1, msLocation = "/usr/local/bin/ms", 
    compareLocation = "comparecladespipe.pl", assign = "assign.txt", 
    observed = "observed.txt", unresolvedTest = TRUE, print.ms.string = FALSE, 
    print.results = FALSE, debug = FALSE, method = "nlminb", 
    itnmax = NULL, return.all = FALSE, maxtime = 0, maxeval = 0, 
    ...) 
{
    modelID <- returnModel(p, migrationArrayMap)
    if (print.results) {
        resultVector <- c(modelID, p)
        names(resultVector) <- c("migrationArryIndividualID", 
            "collapseMatrix.number", "n0multiplierMap.number", 
            "migrationArray.number")
        print(resultVector)
    }
    if (is.na(modelID)) {
        return(badAIC)
    }
    else {
        paramNames <- msIndividualParameters(migrationArray[[modelID]])
        startingVals <- log(c(rlnorm(sum(grepl("collapse", paramNames)), 
            1, 1), rlnorm(sum(grepl("n0multiplier", paramNames)), 
            1, 1), rbeta(sum(grepl("migration", paramNames)), 
            shape1 = 1, shape2 = 3)))
        if (debug) {
            print(startingVals)
        }
        searchResults <- nloptr(x0 = startingVals, eval_f = returnAIC, 
            opts = list(maxeval = itnmax, algorithm = "NLOPT_LN_SBPLX", 
                print_level = 1, maxtime = maxtime, maxeval = maxeval), 
            migrationIndividual = migrationArray[[modelID]], 
            popVector = popVector, badAIC = badAIC, maxParameterValue = maxParameterValue, 
            nTrees = nTrees, msLocation = msLocation, compareLocation = compareLocation, 
            assign = assign, observed = observed, unresolvedTest = unresolvedTest, 
            print.ms.string = print.ms.string, print.results = print.results, 
            debug = debug)
        ifelse(return.all, return(searchResults), return(searchResults$objective))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
