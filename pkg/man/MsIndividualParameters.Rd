\name{MsIndividualParameters}
\alias{MsIndividualParameters}
\title{
MsIndividualParameters
}
\description{
Show the distribution of parameters of a certain model for a given \link{migrationArray}. The considered parameters are: coalescent events (collapse), population size (n0multiplier) or, different migration parameters (migration).
}
\usage{
MsIndividualParameters(migrationIndividual)
}
\author{
Brian O'Meara, Bryan Carstens, Nathan Jackson, Ariadna Morales-Garcia

Maintainer: Brian O'Meara <bomeara@utk.edu>
}
\seealso{
\link{migrationArray}
}
\examples{
	
## In this example model number 54 is shown when there are 3 populations and 3 free
##parameters were specified (maxK=3)  

#load("/path/pkg/data/MigrationArray_3pop_3K.rda")#msIndividualParameters(migrationArray[[54]])
}

