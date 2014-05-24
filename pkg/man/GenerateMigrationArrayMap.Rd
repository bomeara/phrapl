\name{GenerateMigrationArrayMap}
\alias{GenerateMigrationArrayMap}
\title{
GenerateMigrationArrayMap
}
\description{

This function was meant to declare the matrix of model space search for a given \link{migrationArray}. It will generate a table which contains of all the parameter switches for each model.
}
\usage{
GenerateMigrationArrayMap(migrationArray)
}
\arguments{
  \item{migrationArray}{
Contains a matrix of models to be explored.
}
}
\author{
  Brian O'Meara, Bryan Carstens, Nathan Jackson, Ariadna Morales-Garcia

  Maintainer: Brian O'Meara <bomeara@utk.edu>
}
\references{
O'Meara, B.C., N. Jackson, A. Morales-Garcia, and B. Carstens (2014) Phrapl in prep.
}
\seealso{
\link{migrationArray}
}
\examples{
#load("/path/pkg/data/MigrationArray_3pop_3K.rda")
#GenerateMigrationArrayMap(migrationArray)
}


