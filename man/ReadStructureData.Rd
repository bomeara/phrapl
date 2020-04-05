\name{ReadStructureData}
\alias{ReadStructureData}
\title{
Read Structure data
}
\description{
This function reads a Structure file for SNPs. It assumes you have states 1 and
2 and that missing data is encoded as -9.

It returns a list with three objects:

$snps the data to pass into ConvertStructureDataToTrees

$sample.names the first column of sample names

$populations the second column in the original data, assumed (but not required)
to be population names.
}
\usage{
ReadStructureData(file, pairs)
}
\arguments{
  \item{file}{
  A file with the structure data
  }
  \item{pairs}{
  If 0, take all columns and rows. If 1, assume individuals are repeated on rows
  and so you should sample one row at random per site. If 2, assume individuals are
  represented as duplicated columns, so choose samples at random.
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
