\name{kCollapseMatrix}
\alias{kCollapseMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
kCollapseMatrix
}
\description{
Show the number of possible free parameters for times of populations merging.
}
\usage{
kCollapseMatrix(collapseMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collapseMatrix}{
%%     ~~Describe \code{collapseMatrix} here~~
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (collapseMatrix) 
{
    maxVector <- colMax(collapseMatrix)
    return(length(which(maxVector > 0)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
