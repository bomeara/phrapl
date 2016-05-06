\name{CalculateModelAverages}
\alias{CalculateModelAverages}
\title{
Calculate Model Averages
}
\description{
This function outputs model averaged parameter estimates for a set of models.
}
\usage{
CalculateModelAverages(totalData, averageAcrossAllModels = TRUE, parmStartCol = 9)
}
\arguments{
  \item{totalData}{
Path to input file (output of \link{ConcatenateResults})
}
  \item{averageAcrossAllModels}{
If TRUE, each parameter is model averaged across all models in the dataset, even if that model
does not include the given parameter. In these cases where a parameter is absent (i.e., the value 
is NA) this method simply assumes that the value for that parameter is zero. If FALSE, each 
parameter is model averaged by only considering those models that include the relevant parameter.

The former method may be most appropriate for migration, as models that exclude this parameter 
are effectively setting migration to zero. One should be careful interpreting model averaged 
coalescence times (t), however. Keep in mind that when a model excluding a particular coalescent 
event recieves high support, this can be subject to different interpretations. For example, this
could signal that these two populations are so similar that coalescence time between them is
effectively zero. However, it could also signal that these two populations are very distinct, 
but that they coalesce with other populations prior to coalescing with each other.
}
  \item{parmStartCol}{
This gives the column number in totalData in which the parameter values begin. If using default 
PHRAPL output, this should be column 9.
}
}
\seealso{
\link{ConcatenateResults}
}
\examples{
#totalData<-load("/path/totalData.txt")
#modelAverages<-CalculateModelAverages(totalData)
}

