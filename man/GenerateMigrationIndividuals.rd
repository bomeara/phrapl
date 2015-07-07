\name{GenerateMigrationIndividuals}
\alias{GenerateMigrationIndividuals}
\alias{GenerateMigrationIndividual}
\title{
Generate a migrationArray
}
\description{
This function produces an array of all possible demographic models given a number of
populations (specified by popVector), a maximum number of free parameters (maxK), and 
other optional filtering criteria.
}
\usage{
GenerateMigrationIndividuals(popVector, maxK = SetMaxK(popVector), maxMigrationK = 2, 
    maxN0K = 1, forceSymmetricalMigration = TRUE, forceTree = FALSE, 
    verbose = FALSE, parallelRep = NULL) 
}
\arguments{
  \item{popVector}{
A vector that gives the number of samples for each population in a dataset. The model set generated
by this function is based on the number of populations under consideration (indicated by length(popVector)).
}
  \item{maxK}{
The maximum number of free parameters to be incorporated into the model set. The default is calculated based
on the number of samples in popVector.
}
  \item{maxMigrationK}{
The maximum number of migration rate parameters to be incorporated into the model set.
}
  \item{maxN0K}{
The maximum number of n0 parameters to be incorporated into the model set.
}
  \item{forceSymmetricalMigration}{
If TRUE, migration rate indexes between populations are forced to be equal. 
}
  \item{forceTree}{
If TRUE, only models with fully resolved topologies are included in the model set.
}
  \item{verbose}{
If TRUE, prints out status updates
}
  \item{parallelRep}{
If running the function through loop in parallel on a cluster, this argument is assigned 
a replicate number.
}
}
\author{
Brian O'Meara
}
\note{
For more information, please see the user manual.
}
\keyword{ ~migrationIndividual }
