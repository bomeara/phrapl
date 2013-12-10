\name{setMaxK}
\alias{setMaxK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
setMaxK
}
\description{
Maximum value of free parameters to be used in the searching process
}
\usage{
setMaxK(popVector, nTrees = 1, taxaPerParam = 5, force = FALSE)
}
\arguments{
  \item{popVector}{
%%     ~~Describe \code{popVector} here~~
}
  \item{nTrees}{
%%     ~~Describe \code{nTrees} here~~
}
  \item{taxaPerParam}{
%%     ~~Describe \code{taxaPerParam} here~~
}
  \item{force}{
%%     ~~Describe \code{force} here~~
}
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (popVector, nTrees = 1, taxaPerParam = 5, force = FALSE) 
{
    if (exists("maxK") && !force) {
        print(paste("note: setMaxK called, but maxK already exists, so using the existing one of ", 
            maxK, ". To change this behavior, give setMaxK the option force=TRUE", 
            sep = ""))
        return(maxK)
    }
    else {
        maxK <- max(1, floor(sum(popVector)/(nTrees * taxaPerParam)))
        return(maxK)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
