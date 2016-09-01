\name{CheckMS}
\alias{CheckMS}
\alias{CheckMs}
\title{
A check of whether the pre-compiled version of ms is compatible with your operating system
}
\description{
We have packaged the raw code for the program ms (Hudson 2002) with PHRAPL, which PHRAPL uses to simulate gene trees. 
One can compile this code upon initally downloading PHRAPL using the \code{CompileMS} function. However, we have also 
included an ms version that has been pre-complied for a Mac. The purpose of \code{CheckMS} is to check whether 
this pre-compiled version is compatible with your operating system. If running this function returns a message that 
ms is working, then there is no need to compile ms prior to running PHRAPL. If the function returns an error or
a message that there is something wrong, then you will need to compile ms using the \code{CompileMS} function. Once
you have compiled ms, this check can again be run to ensure that the compilation was sucessful.
}
\usage{
CheckMS()
}
\references{
Hudson, R.R. 2002. Generating samples under a Wright-Fisher neutral model. \emph{Bioinformatics}. 18:337-338.
\url{http://home.uchicago.edu/~rhudson1/source.html}
}
\author{
Nathan Jackson
}