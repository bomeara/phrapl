\name{MsIndividualParameters}
\alias{MsIndividualParameters}
\title{
MsIndividualParameters
}
\description{
Show the distribution of parameters of a certain model for a given
\link{migrationArray}. The considered parameters are: coalescent events 
(collapse), population size (n0multiplier) or, different migration parameters 
(migration).
}
\usage{
MsIndividualParameters(migrationIndividual)
}
\arguments{
  \item{migrationIndividual}{
This is the chosen model. If you want model 17 from migrationArray, you could 
pass this in as migrationIndividual=migrationArray[[17]]
}
}
\author{
Brian O'Meara, Bryan Carstens, Nathan Jackson, Ariadna Morales-Garcia

Maintainer: Brian O'Meara <bomeara@utk.edu>
}
\seealso{
\link{migrationArray}
}
\examples{
	
## In this example model number 54 is shown when there are 3 populations and 3 
##free parameters were specified (maxK=3)  

data(MigrationArray_4pop_3K)
msIndividualParameters(migrationArray[[54]])
}

