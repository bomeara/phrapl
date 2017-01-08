\name{ConvertStructureDataToTrees}
\alias{ConvertStructureDataToTrees}
\title{
Use Structure data
}
\description{
This function takes the snps output of ReadStructureData() and creates a set of
trees you can use in GridSearch or elsewhere with the doSNP=TRUE argument set.
It returns a multiphylo object. Missing data (-9 in the original file) are removed.
}
\usage{
ConvertStructureDataToTrees(snps)
}
\arguments{
  \item{snps}{
The output of ReadStructureData, the snps object
}
}
\author{
  Brian O'Meara, Bryan Carstens, Nathan Jackson, Ariadna Morales-Garcia

  Maintainer: Brian O'Meara <bomeara@utk.edu>
}
\references{
  O'Meara, B.C., N. Jackson, A. Morales-Garcia, and B. Carstens (2014)
  Phrapl in prep.
}

\examples{

#snptrees <- ConvertStructureDataToTrees(ReadStructureData(file="mydata")$snps)

}
