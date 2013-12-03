\name{returnAIC}
\alias{returnAIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
returnAIC
}
\description{
Show the AIC score for the specified model.
}
\usage{
returnAIC(par, popVector, migrationIndividual, nTrees = 1, msLocation = "/usr/local/bin/ms", compareLocation = "comparecladespipe.pl", assign = "assign.txt", observed = "observed.txt", unresolvedTest = TRUE, print.results = FALSE, print.ms.string = FALSE, debug = FALSE, badAIC = 1e+14, maxParameterValue = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
%%     ~~Describe \code{par} here~~
}
  \item{popVector}{
%%     ~~Describe \code{popVector} here~~
}
  \item{migrationIndividual}{
%%     ~~Describe \code{migrationIndividual} here~~
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
  \item{print.results}{
%%     ~~Describe \code{print.results} here~~
}
  \item{print.ms.string}{
%%     ~~Describe \code{print.ms.string} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
}
  \item{badAIC}{
%%     ~~Describe \code{badAIC} here~~
}
  \item{maxParameterValue}{
%%     ~~Describe \code{maxParameterValue} here~~
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
function (par, popVector, migrationIndividual, nTrees = 1, msLocation = "/usr/local/bin/ms", 
    compareLocation = "comparecladespipe.pl", assign = "assign.txt", 
    observed = "observed.txt", unresolvedTest = TRUE, print.results = FALSE, 
    print.ms.string = FALSE, debug = FALSE, badAIC = 1e+14, maxParameterValue = 100) 
{
    parameterVector <- exp(par)
    names(parameterVector) <- msIndividualParameters(migrationIndividual)
    if (max(parameterVector) > maxParameterValue) {
        return(badAIC)
    }
    if (debug) {
        print(parameterVector)
    }
    print("wd is ")
    print(getwd())
    likelihoodVector <- pipeMS(popVector = popVector, migrationIndividual = migrationIndividual, 
        parameterVector = parameterVector, nTrees = nTrees, msLocation = msLocation, 
        compareLocation = compareLocation, assign = assign, observed = observed, 
        unresolvedTest = unresolvedTest, print.ms.string = print.ms.string, 
        debug = debug)
    lnLValue <- convertOutputVectorToLikelihood(likelihoodVector, 
        nTrees = nTrees, probOfMissing = 1/howmanytrees(sum(popVector)))
    AICValue <- 2 * (-lnLValue + kAll(migrationIndividual))
    if (print.results) {
        resultsVector <- c(AICValue, lnLValue, exp(par))
        names(resultsVector) <- c("AIC", "lnL", msIndividualParameters(migrationIndividual))
        print(resultsVector)
    }
    return(AICValue)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
