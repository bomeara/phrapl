\name{GenerateMigrationIndividualsOneAtATime}
\alias{GenerateMigrationIndividualsOneAtATime}
\alias{GenerateMigrationIndividualOneAtATime}
\title{
GenerateMigrationIndividualsOneAtATime
}
\description{
The purpose of this function is create a single a priori migrationIndividual
for a given collapse, n0multiplier, growth, and migration history.
The four inputs for this function are collapseList, n0multiplierMap, 
growthMap, and migrationList. 

CollapseList is the only set of parameters that must be specified. This gives a
list of collapse history vectors, one vector for each 
coalescent event in the tree. So, collapseList = list(c(1,1,0),c(1,NA,1)) means 
that there are two coalescent events: in the first event, population 1 and 2 
coalesce while population 3 does not; in the second event, ancestral population 
1-2 coalesces with population 3. 

The remaining three parameters may be specified or not specified. If not specified,
(i.e., set to NULL),then null matrices will be automatically constructed in which 
all n0multipliers are set to one, and all growth and migration parameters are set to zero.

If specifying n0multiplierList and/or growthList, the format for these is the same as for
collapseList, and the available parameters for these must match the splitting history 
depicted in the collapseList.

MigrationList is a list of migration matrices. There will be one migration matrix 
for the present generation, plus any historical matrices that apply 
(there will be as many matrices as there are collapse events).
So, in the three population scenario, if there is symmetrical migration between 
populations 1 and 2 and no historical migration, migrationList will be:

migrationList<-list(\cr
 	t(array(c(\cr
	NA, 1, 0,\cr
	1, NA, 0,\cr
	0, 0, NA),\cr
 	dim=c(3,3))),
 
 	t(array(c(\cr
	NA, NA, 0,\cr
	NA, NA, NA,\cr
	0,  NA, NA),\cr
 	dim=c(3,3))))

Note that in R, arrays are constructed by reading in values from column 1 first then from 
column 2, etc. However, it is more intuitive to construct migration matrices by rows (i.e., 
first listing values for row 1, then row 2, etc). Thus, I think it is easier to type in the 
arrays by rows (as done above), and then transpose them (using "t"). Also, spacing and hard 
returns can be used to visualize these values in the form of matrices.
}
\usage{
GenerateMigrationIndividualsOneAtATime(collapseList=collapseList,n0multiplierList=n0multiplierList,
	growthList=growthList,migrationList=migrationList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collapseList}{
List of collapse histories in the model
}
  \item{n0multiplierList}{
List of population size scalars in the model
}
  \item{growthList}{
List of growth rate scalars in the model
}
  \item{migrationList}{
List of migration rate matrices in the model
}
}
\author{
  Nathan Jackson and Brian O'Meara
}
\references{
  O'Meara, B.C., N. Jackson, A. Morales-Garcia, and B. Carstens (2014) 
  Phrapl in prep.
}
\examples{

# ##Four population Example (with coalescence, growth, and migration parameters)
# 
# #Create list of collapse histories
# collapse_1<-c(0,0,1,1)
# collapse_2<-c(1,1,0,NA)
# collapse_3<-c(1,NA,1,NA)
# collapseList<-list(collapse_1,collapse_2,collapse_3)
# 
# #Create list of growth histories
# growth_1<-c(0,0,1,1)
# growth_2<-c(0,0,0,NA)
# growth_3<-c(0,NA,0,NA)
# growthList<-list(growth_1,growth_2,growth_3)
# 
# #Create list of migration matrices
# migration_1<-t(array(c(
# 	NA, 0, 0, 0,
# 	0, NA, 0, 0,
# 	0, 0, NA, 1,
# 	0, 0, 1, NA),
# 	dim=c(4,4)))
# 
# migration_2<-t(array(c(
# 	NA, 0, 0, NA,
# 	1, NA, 0, NA,
# 	0, 0, NA, NA,
# 	NA,NA,NA, NA),
# 	dim=c(4,4)))
# 
# migration_3<-t(array(c(
# 	NA, NA, 0, NA,
# 	NA, NA, NA,NA,
# 	0,  NA, NA,NA,
# 	NA, NA, NA,NA),
# 	dim=c(4,4)))
# 
# migrationList<-list(migration_1,migration_2,migration_3)
# 
# #Run function
# migrationIndividual<-GenerateMigrationIndividualsOneAtATime(collapseList=collapseList,growthList=growthList,
# 	migrationList=migrationList)



# ##Three population Example (with coalescence and migration parameters)
# 
# #Create list of collapse histories
# collapse_1<-c(1,1,0)
# collapse_2<-c(1,NA,1)
# collapseList<-list(collapse_1,collapse_2)
# 
# #Create list of migration matrices
# migration_1<-t(array(c(
# 	NA, 2, 1,
# 	2, NA, 1,
# 	1, 1, NA),
# 	dim=c(3,3)))
# 
# migration_2<-t(array(c(
# 	NA, NA, 3,
# 	NA, NA, NA,
# 	3,  NA, NA),
# 	dim=c(3,3)))
# 
# migrationList<-list(migration_1,migration_2)
# 
# #Run function
# migrationIndividual<-GenerateMigrationIndividualsOneAtATime(collapseList=collapseList,
# 	migrationList=migrationList)

}
