\name{CreateAssignment}
\alias{CreateAssignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CreateAssignment
}
\description{
Define the number of free parameters and assign individuals to populations according to the predefined file of samples to populations.
}
\usage{
CreateAssignment(popVector, file = "assign.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popVector}{
%%     ~~Describe \code{popVector} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
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
function (popVector, file = "assign.txt") 
{
    alphabet <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ", "")[[1]]
    assignFrame <- data.frame()
    indivTotal <- 0
    for (popIndex in 1:length(popVector)) {
        popLabel <- alphabet[popIndex]
        for (indivIndex in 1:popVector[popIndex]) {
            indivTotal <- indivTotal + 1
            if (indivTotal == 1) {
                assignFrame <- data.frame(popLabel, indivTotal)
            }
            else {
                assignFrame <- rbind(assignFrame, data.frame(popLabel, 
                  indivTotal))
            }
        }
    }
    write.table(assignFrame, file = file, quote = FALSE, sep = "\t", 
        row.names = FALSE, col.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
