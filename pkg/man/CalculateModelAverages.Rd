\name{CalculateModelAverages}
\alias{CalculateModelAverages}
\title{
CalculateModelAverages
}
\description{
This function outputs model averaged parameter estimates
}
\usage{
CalculateModelAverages(totalData = totalData, parmStartCol = 10)
}
\arguments{
  \item{totalData}{
Path to input file (output of \link{ConcatenateResults})
}
  \item{parmStartCol}{
The column number in which the parameter values begin (col 9)
}
}
\seealso{
\link{ConcatenateResults}
}
\examples{
#totalData<-load("/path/totalData.txt")
#modelAverages<-CalculateModelAverages(totalData, parmStartCol=9)
}

