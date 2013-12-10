\name{returnModel}
\alias{returnModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
returnModel
}
\description{
Show a list of the models that match better to the data
}
\usage{
returnModel(p, migrationArrayMap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{migrationArrayMap}{
%%     ~~Describe \code{migrationArrayMap} here~~
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
function (p, migrationArrayMap) 
{
    prunedResults <- subset(migrationArrayMap, migrationArrayMap$collapseMatrix.number == 
        p[1])
    prunedResults <- subset(prunedResults, prunedResults$n0multiplierMap.number == 
        p[2])
    prunedResults <- subset(prunedResults, prunedResults$migrationArray.number == 
        p[3])
    if (dim(prunedResults)[1] == 1) {
        return(prunedResults$model)
    }
    else {
        warning(paste("All these models match", prunedResults$model, 
            "taking first one"))
        print(paste("All these models match", prunedResults$model, 
            "taking first one"))
        return(prunedResults$model[1])
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
