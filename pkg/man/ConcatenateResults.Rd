\name{ConcatenateResults}
\alias{ConcatenateResults}
\alias{concatenateResults}
\alias{concatenate.Results}
\title{
Cull results from a grid or heuristic search into a table
}
\description{
This function takes output from a grid or heuristic search carried out across a set of models and constructs a 
table that includes AIC values and weights, parameter estimates, and specifications for each model. The function 
can either import one or more R object files (with the extention .rda) located in the specified rdaFilePath, in 
which phrapl output is stored, or can can be passed a vector of file names to be imported.
}
\usage{
ConcatenateResults(rdaFilesPath="./",rdaFiles=NULL,outFile=NULL,rmNaParameters=TRUE)
}
\arguments{
  \item{rdaFilesPath}{
A directory in which R object files (which must have the extention .rda) are located. 
}
  \item{rdaFiles}{
As an alternative to using rdaFilePath, a vector of file names can be passed to the function.
}
  \item{outFile}{
If specified, gives the file name to which the output table should be printed.
}
  \item{rmNaParameters}{
If true, this removes all parameter columns from the table that only include NAs.
}
\author{
Nathan Jackson
}
\note{
For more information, please see the user manual.
}
\keyword{ ~concatenate }
