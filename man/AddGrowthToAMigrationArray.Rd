\name{AddGrowthToAMigrationArray}
\alias{AddGrowthToAMigrationArray}
\title{
Add growth paramter matrices (growthMaps) to a list of models (a migrationArray) that
lack growth matrices
}
\description{
Previous versions of PHRAPL produced model lists (migrationArrays) that lack growth
matrices (growthMaps). However, growthMaps must now be present to analyze models in PHRAPL, 
even if growth parameters are not incorporated in the model. Thus, to facilitate the analysis 
of old migrationArrays, this function takes as input a  migrationArray lacking growthMaps
and adds an empty growthMap to each model (i.e., a growthMap filled with zeros) within the 
migrationArray.
}
\usage{
AddGrowthToAMigrationArray(migrationArray)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{migrationArray}{
A list of models to which NULL growth matrices should be added (i.e., matrices positing zero
growth). A single model can also be specified.
}
}
\author{
  Nathan Jackson
}
\references{
  O'Meara, B.C., N. Jackson, A. Morales-Garcia, and B. Carstens (2014) 
  Phrapl in prep.
}
