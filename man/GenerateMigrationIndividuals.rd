\name{GenerateMigrationIndividuals}
\alias{GenerateMigrationIndividuals}
\title{
Generate a list of models (a migrationArray) to analyze using PHRAPL
}
\description{
This function produces a list of all possible demographic models given a number of populations
(specified by popVector), number of free parameters (specified by maxK), and other optional
filtering criteria. The list of models is call a migrationArray (a name which, confusingly,
is also given to migration matrices within a model) and each model is referred to as a
migrationIndividual.

Each migrationIndividual is composed of four major demographic components:
1. The coalescence history (the collapseMatrix, in which rows correspond to populations and columns
	correspond to temporal events)
2. Population size scalar parameters (the n0multiplierMap, with the same dimensions as CollapseMatrix)
3. Growth (alpha) parameters (the growthMap, also with the same dimensions as CollapseMatrix)
4. Migration rate matrices (the migrationArray, where the number of matrices is equal to the number
	of temporal events specified by the collapseMatrix.

One can specify the total maximum number of free parameters (maxK), and also specify the maximum
number of free parameters for a given parameter type (e.g., maxN0K, maxGrowthK, and maxMigrationK).
Note that maxGrowthK and maxMigrationK can be set to zero; however the lowest possible value for
N0K is one (in which case all population sizes are the same).

Specific demographic components can also be fixed, such that the generated migrationArray only varies
subsets of parameters. Fixed components are specified by collapseList, n0multiplierList, growthList,
and migrationList. Note that a collapseList must always be specified in order to specify any other
demographic component. The format for specifying collapseList, n0multiplierList, and growthList are
the same: as a list of vectors. A vector contains parameter values for each population (e.g.,
becoming the rows in collapseMatrix), and each vector in the list represents a different temporal
event (e.g., becoming the columns in collapseMatrix).
For example, collapseList = list(c(1,1,0),c(2,NA,2)) means that there are two coalescent events:
in the first event, population 1 and 2 coalesce while population 3 does not; in the second event,
ancestral population 1-2 coalesces with population 3.

MigrationList differs in that a list of matrices, rather than vectors, must be specified. There will
be one migration matrix for the present generation, plus any historical matrices that apply (there
will be as many matrices as there are collapse events). So, in the three population scenario, if there
is symmetrical migration between populations 1 and 2 and no historical migration, migrationList will be:

 migrationList<-list(\cr
 t(array(c(\cr
 NA, 1, 0,\cr
 1, NA, 0,\cr
 0, 0, NA),\cr
 dim=c(3,3))),\cr

 t(array(c(\cr
 NA, NA, 0,\cr
 NA, NA, NA,\cr
 0,  NA, NA),\cr
 dim=c(3,3))))\cr

Note that in R, arrays are constructed by reading in values from column 1 first then from column 2, etc.
However, it is more intuitive to construct migration matrices by rows (i.e., first listing values for
row 1, then row 2, etc). Thus, in the example above, arrays are entered as rows, and then transposed
(using "t"). Also, spacing and hard returns can be used to visualize these values in the form of matrices.

Other methods of model filtering (using forceSymmetricalMigration or forceTree) can also be implemented.
}
\usage{
GenerateMigrationIndividuals(popVector,maxK=SetMaxK(popVector),maxN0K=1,maxGrowthK=0,
	maxMigrationK=1,collapseList=NULL,n0multiplierList=NULL,growthList=NULL,migrationList=NULL,
	forceSymmetricalMigration=TRUE,forceTree=FALSE,verbose=FALSE,parallelRep=NULL)
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
  \item{maxN0K}{
The maximum number of n0 parameters to be incorporated into the model set.
}
  \item{maxGrowthK}{
The maximum number of growth parameters to be incorporated into the model set.
}
  \item{maxMigrationK}{
The maximum number of migration rate parameters to be incorporated into the model set.
}
  \item{collapseList}{
A particular collapse history to be specified
}
  \item{n0multiplierList}{
A particular n0muliplier history to be specified
}
  \item{growthList}{
A particular population growth history to be specified
}
  \item{migrationList}{
A particular migration history to be specified
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
Brian O'Meara and Nathan Jackson
}
\note{
For more information, please see the user manual.
}
\keyword{ ~migrationIndividual }
\examples{

# #Assuming a given tree and migration scenario, create a model set that tests all
# #possible distributions of a single population growth parameter across three populations.
#
# popVector<-c(3,3,3)
# maxK=5
# maxGrowthK=2
# forceTree=TRUE
#
# #Fix a tree
# collapse_1<-c(1,1,0)
# collapse_2<-c(2,NA,2)
# collapseList<-list(collapse_1,collapse_2)
#
# #Fix migration
# migration_1<-t(array(c(
# 	NA, 1, 1,
# 	1, NA, 1,
# 	1, 1, NA),
# dim=c(3,3)))
#
# migration_2<-t(array(c(
# 	NA, NA, 2,
# 	NA, NA, NA,
# 	2,  NA, NA),
# dim=c(3,3)))
# migrationList<-list(migration_1,migration_2)
#
# migrationArray<-GenerateMigrationIndividuals(popVector=popVector,maxK=maxK,maxGrowthK=maxGrowthK,
# 	collapseList=collapseList,migrationList=migrationList,forceTree=forceTree)

}
