\name{ConcatenateResults}
\alias{ConcatenateResults}
\title{
ConcatenateResults
}
\description{
This function builds a table with details for each model analyzed. It give AICs, model weights, details of each model, and parameter estimates and indexes. This is the major output of phrapl.
}
\usage{
ConcatenateResults(rdaFilesPath = "./", rdaFiles = NULL, outFile = NULL, addAICweights = TRUE, rmNaParameters = TRUE, addTime.elapsed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdaFilesPath}{
Path to *.rda files which are outputs (not automatically saved!) of \link{GridSearch} (default is the current dir).

}
  \item{rdaFiles}{
A vector of rda file names (located in the directory specified above) can be given.
}
  \item{outFile}{
Outfile name.
}
  \item{addAICweights}{
Option to printed out in outFile.
}
  \item{rmNaParameters}{
Option to printed out in outFile.
}
  \item{addTime.elapsed}{
Option to printed out in outFile.
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
\link{GridSearch}, 
\link{CalculateModelAverages}.
}
\examples{

#totalData<-ConcatenateResults(
#rdaFilesPath="/path_to_rdaFiles/",
#rdaFiles=NULL,
#outFile=NULL,
#addAICweights=TRUE,
#rmNaParameters=TRUE,
#addTime.elapsed=FALSE)

##Save results to a .txt file
#write(t(totalData), file = "totalData.txt", ncolumns = 40, append = FALSE, sep = "\t")	
		
}

