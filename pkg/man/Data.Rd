\name{phrapl-data}
\alias{MigrationArray_3pop_3K}
\alias{MigrationArray_4pop_3K}
\alias{migrationArray_3pop_3K}
\alias{migrationArray_4pop_3K}
\title{
Cached migrationArrays
}
\description{
These are precomputed migrationArrays and migrationArrayMaps (soon). 

To access one with 3 populations and no more than 3 free parameters, use
data(MigrationArray_3pop_3K)

To access one with 4 populations and no more than 3 free parameters, use
data(MigrationArray_4pop_3K)


Two objects are included in each 
1. migrationArray, which lists specifications for  models to be analyzed 
2. migrationArrayMap, which lists parameter indexes for each model included 
in migrationArray (not yet available)

}
\examples{
  data(MigrationArray_3pop_3K)
  length(migrationArray)
}
\keyword{example}
\keyword{TestData}
