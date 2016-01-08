\name{AddEventToMigrationArray}
\alias{AddEventToMigrationArray}
\alias{AddEventToAMigrationArray}
\title{
Add an Event To a MigrationArray
}
\description{
This function integrates a non-coalescence demographic event within a model
or set of models. This can be useful if one wants to posit shifts in a parameter that
do not correspond to a splitting event (e.g., one would like migration to only occur 
for a given length of time after a split, but then to cease).
}
\usage{
AddEventToMigrationArray<-function(migrationArray,eventTime,n0multiplierVec=NULL,growthVec=NULL,migrationMat=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{migrationArray}{
A list of models (or a single model) to which an event should be added
}
  \item{eventTime}{
The relative time period at which the new event should be inserted within the collapseMatrix
}
  \item{n0multiplierVec}{
A vector or single value specifying the n0multiplier parameter indexes to be invoked during 
the new time interval. If \code{NULL}, no new set of parameter indexes are invoked at this time.
}
  \item{growthVec}{
A vector or single value specifying the growth parameter indexes to be invoked during 
the new time interval. If \code{NULL}, no new set of parameter indexes are invoked at this time.
}
  \item{migrationMat}{
A matrix, vector, or single value specifying the migration parameter indexes to be invoked during 
the new time interval. If \code{NULL}, no new set of parameter indexes are invoked at this time.
}
}
\details{
To use this function, one must specify a model (migrationIndividual) or set of models
(migrationArray). If a set of models is specified, these models must all contain the 
same number of populations and the same number of collapse events (i.e., the collapseMatrix
compenent within each migrationIndividual must have the same dimensions). 

The relative timing of the desired new event must be specified as a single number using 
eventTime. An eventTime of 1 will place a new event (i.e., column) in the first column 
position within the collapseMatrix (and the other columns will be shifted to the right). 
An eventTime of 2 will place the new column in the second position, etc. The eventTime
cannot exceed the number of events (i.e., columns) within the original collapseMatrix.
The added column will consist entierly of NAs, indicating that no population coalescence 
can occur at this time.

Finally, one can specify a new set of n0multiplier, growth, and/or migration parameter
indexes to be invoked at the new specified time period using n0muliplierVec, growthVec, or
migrationMat, respectively. When the default value for these is used (which is \code{NULL}), 
n0multiplier, growth, and migration matrices within each model are automatically expanded 
by simply copying over parameter indexes from the adjacent time period to the new time
period (i.e., no change is invoked). For n0multiplier and growth, a new column is added; 
for migration, a new matrix is added. 

However, one can also specify the parameter indexes to be used at this new time, allowing
a shift in that parameter, without a correponding coalescence event. These indexes can be
specified in one of two ways:

First, one can simply specify a single parameter index which will be implemented for all 
relevant populations at that time period. For example, if the following collapseMatrix:

\code{$collapseMatrix\cr
	  [,1] [,2]\cr
[1,]    1    2\cr
[2,]    1   NA\cr
[3,]    0    2
}

is expanded to include a shift from non-zero migration among all populations to zero migration 
at some point prior to the coallesce of populations 1 and 2 (meaning that migration has 
occurred recently--upon secondary contact--but did not accompany initial divergence), we 
should set \code{eventTime = 1} and \code{migrationMat = 1}. Assuming that the original model included 
zero migration across all time periods, this will produce the following collpaseMatrix:

\code{$collapseMatrix\cr
     [,1] [,2] [,3]\cr
[1,]   NA	 1    2\cr
[2,]   NA	 1   NA\cr
[3,]   NA	 0    2
}

And the migrationArray will be expanded from two matrices:

\code{$migrationArray\cr
, , 1\cr
\cr
     [,1] [,2] [,3]\cr
[1,]   NA    0    0\cr
[2,]    0   NA    0\cr
[3,]    0    0   NA\cr
\cr
, , 2\cr
\cr
     [,1] [,2] [,3]\cr
[1,]   NA   NA    0\cr
[2,]   NA   NA   NA\cr
[3,]    0   NA   NA
}


to three matrices:

\code{$migrationArray\cr
, , 1\cr
\cr
     [,1] [,2] [,3]\cr
[1,]   NA    1    1\cr
[2,]    1   NA    1\cr
[3,]    1    1   NA\cr
\cr
, , 2\cr
\cr
     [,1] [,2] [,3]\cr
[1,]   NA    0    0\cr
[2,]    0   NA    0\cr
[3,]    0    0   NA\cr
\cr
, , 3\cr
\cr
     [,1] [,2] [,3]\cr
[1,]   NA   NA    0\cr
[2,]   NA   NA   NA\cr
[3,]    0   NA   NA
}

where the new migration matrix was added in the first time period position and was filled with 
zeros, as per the value specified.

This approach works similarly for n0multiplier and growth. For example, if one would like to model 
population growth during the initial stages of divergence, followed by population size stability,
one should set \code{eventTime = 1} and \code{growthVec = 0} to change this history:

\code{$collapseMatrix\cr
     [,1] [,2]\cr
[1,]    1    2\cr
[2,]    1   NA\cr
[3,]    0    2\cr
\cr
$growthMap\cr
     [,1] [,2]\cr
[1,]    1    0\cr
[2,]    1   NA\cr
[3,]    0    0
}


into this history:

\code{$collapseMatrix\cr
      [,1] [,2] [,3]\cr
[1,]   NA	 1    2\cr
[2,]   NA	 1   NA\cr
[3,]   NA	 0    2\cr
\cr
$growthMap\cr
      [,1] [,2] [,3]\cr
[1,]    0	 1    0\cr
[2,]    0	 1   NA\cr
[3,]    0	 0    0
}

If one wishes to specify more complex histories, in which different populations have different
indexes, one can alternatively specify the entire new vector (for n0multiplier and growth) 
or matrix (for migration) of parameter indexes to be inserted. For example, in the previous example, 
if one would rather model population size stability upon divergence of populations 1 and 2 followed 
by recent growth of population 1, but not of populations 2 or 3, one can set \code{eventTime = 1} and 
\code{growthVec = c(1,0,0)} to change the growth history from this:

\code{$growthMap\cr
     [,1] [,2]\cr
[1,]    0    0\cr
[2,]    0   NA\cr
[3,]    0    0
}


to this:

\code{$growthMap\cr
      [,1] [,2] [,3]\cr
[1,]    1	 0    0\cr
[2,]    0	 0   NA\cr
[3,]    0	 0    0
}

One can also specify a migration matrix using migrationMat. For example, if, one would like to 
insert a new matrix in which the migration rate differs among population pairs, one would specify

\code{migrationMat = t(array(c(\cr
	NA, 3, 2,\cr
	1, NA, 1,\cr
	2, 3, NA),\cr
dim=c(3,3)))
}
		
If one would rather, this can also be specified simply as a vector (with or without the NAs). For
example, the above migration matrix could be also specified as either \code{c(NA,1,2,3,NA,3,2,1,NA)} 
or as \code{c(1,2,3,3,2,1)}, where the vector is read first by column, then by row. Note that the 
number of values must match the number expected for that time period. For example, if one would like 
to insert a new migration history at ancestral time period 2, one should remember that there are 
fewer populations at that time period, and thus the number of values specified must be reduced 
accordingly. For example, specifying distinct migration parameters from population 1-2 to 3 versus 
from population 3 to 1-2 would be specified as either

\code{migrationMat = t(array(c(\cr
 	NA, NA, 2,\cr
 	NA, NA, NA,\cr
 	1,  NA, NA),\cr
 dim=c(3,3))) 
}
 			
or as

\code{migrationMat = (1,2)
}

Specifying two distinct population size parameters for ancestral population 1-2 and population 3 at
the second time period would be done by setting \code{eventTime = 2} and \code{n0muliplierVec = c(1,2)}, 
to convert

\code{$n0multiplierMap\cr
     [,1] [,2]\cr
[1,]    1    1\cr
[2,]    1   NA\cr
[3,]    1    1
}


into

\code{$n0multiplierMap\cr
      [,1] [,2] [,3]\cr
[1,]    1	 1    1\cr
[2,]    1	NA   NA\cr
[3,]    1	 2    1
}

Finally, note that although the function can only add a single time period to the collapseMatrix, one can
add multiple non-coalescence time periods by repetedly running the function on the same migrationArray.
}
\note{
When analyzing a migrationArray using \code{\link{GridSearch}} that includes new time periods
produced by \code{AddEventToMigrationArray}, a time value must be specified using the \code{addedEventTime} 
argument. \code{addedEventTime} takes either a single value or vector of values, such that there are as 
many values as added events in the collapseMatrix (sorted from most recent to most ancient).
Time values can be specifed as either a specific time period (in units of 4Ne) or as a relative time
period. When the latter is desired, the \code{addedEventTimeAsScalar} argument in \code{\link{GridSearch}} 
must be set to \code{TRUE} (this is the default), which will cause the values specified by 
\code{addedEventTime} to be treated as scalars, rather than as absolute. Scalar values must be between 0 and 1, 
and will be multiplied by whatever the next collapse time estimate happens to be (which will typically vary across 
models). For example, setting \code{addedEventTimeAsScalar = TRUE} and \code{addedEventTime = 0.1} for the 
following model:

\code{$collapseMatrix\cr
      [,1] [,2] [,3]\cr
[1,]   NA	 1    2\cr
[2,]   NA	 1   NA\cr
[3,]   NA	 0    2\cr
\cr
$growthMap\cr
      [,1] [,2] [,3]\cr
[1,]    0	 1    0\cr
[2,]    0	 0   NA\cr
[3,]    0	 0    0
}

will cause (looking backwards in time) population growth to commence in population 1 after 10\% of the 
branch length leading to coalescence of populations 1 and 2 has been reached.

One final thing to highlight is that the specified timing of these events must fit the modeled chronology of 
coalescent and non-coalescence events in the collapseMatrix. That is, the timing for an added non-coalescence 
event given by \code{addedEventTime} must be later than the collapse event directly preceding the added event
and must be earlier than the collapse event directly following it. For this reason, it is recommended 
that \code{addedEventTimes} be entered as scalars, as this guarantees that the time will always be
earlier than the collapse event that follows (but not that the time will occur after the previous collapse
event). Absolute values can also be used, but remember that in a grid search, different values
are explored, so many parameter combinations in the grid may not fit the specified event times. However, 
PHRAPL filters out all parameter combinations in the parameter grid that do not adhere to the specified
chronology of events in the model. Thus, setting absolute dates may reduce the parameter space explored 
for a model, but this parameter space will always fit the model. 
}
\author{
  Nathan Jackson
}
\references{
  O'Meara, B.C., N. Jackson, A. Morales-Garcia, and B. Carstens (2014) 
  Phrapl in prep.
}
\examples{

# ##First, generate a migrationIndividual that models a specific coalescence history among four populations
# collapse_1<-c(0,0,1,1)
# collapse_2<-c(2,2,0,NA)
# collapse_3<-c(3,NA,3,NA)
# collapseList<-list(collapse_1,collapse_2,collapse_3)
# migrationIndividual<-GenerateMigrationIndividualsOneAtATime(collapseList=collapseList)
# 
# ##Then, add a recent migration (secondary contact) event between sister populations 1 and 2
# migrationArray<-migrationIndividual
# eventTime=1
# migrationMat=c(1,0,0,1,0,0,0,0,0,0,0,0)
# 
# migrationArray<-AddEventToMigrationArray(migrationArray=migrationArray,eventTime=eventTime,
# 	migrationMat=migrationMat)

}
