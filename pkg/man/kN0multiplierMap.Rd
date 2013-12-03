\name{kN0multiplierMap}
\alias{kN0multiplierMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
kN0multiplierMap
}
\description{
Show the value of possible population growth rates.
}
\usage{
kN0multiplierMap(n0multiplierMap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n0multiplierMap}{
%%     ~~Describe \code{n0multiplierMap} here~~
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
function (n0multiplierMap) 
{
    return(max(n0multiplierMap, na.rm = TRUE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
