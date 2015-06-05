\name{ConcatenateResults}
\alias{ConcatenateResults}
\title{
ConcatenateResults
}
\description{
This function builds a table with details for each model analyzed using a grid or heuristic search. 
It give AICs, model weights, details of each model, and parameter estimates and indexes. 
The function can either import one or more R object files (with the extention .rda) located in the 
specified rdaFilePath, in which phrapl output is stored, or can be passed a vector of file names to be 
imported.This is the major output of phrapl.
}
\usage{
ConcatenateResults(rdaFilesPath = "./", rdaFiles = NULL, outFile = NULL, 
addAICweights = TRUE, rmNaParameters = TRUE, addTime.elapsed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdaFilesPath}{
Path to where R object files (which must have the extention .rda) are located which contain
results from \link{GridSearch} (default is the current dir).
}
  \item{rdaFiles}{
As an alternative to reading in all .rda files in the specified directory, a vector of file names (located in the 
directory specified above) can be given.
}
  \item{outFile}{
If specified, gives the file name to which an output table will be printed.
}
  \item{addAICweights}{
If TRUE, will add AICweights to the table.
}
  \item{rmNaParameters}{
If true, this removes all parameter columns from the table that only include NAs.
}
  \item{addTime.elapsed}{
If elapsed time is being tracked, this adds a column for this.
}
}
\details{
If subsets of models from a migrationArray are run in separate analyses, one should save the results as
R object files (.rda) after each analysis (e.g., save(list=ls(), FILENAME)). Once all models are run,
ConcatenateResults can then cull the results from all the files into a master result table.
}
\author{
  Brian O'Meara, Bryan Carstens, Nathan Jackson, Ariadna Morales-Garcia

  Maintainer: Brian O'Meara <bomeara@utk.edu>
}
\references{
 
O'Meara, B.C., N. Jackson, A. Morales-Garcia, and B. Carstens (2014) Phrapl in 
prep.
}
\seealso{
\link{GridSearch}, 
\link{CalculateModelAverages}.
}
\examples{
\dontrun{
totalData<-ConcatenateResults(
rdaFilesPath="/path_to_rdaFiles/",
rdaFiles=NULL,
outFile=NULL,
addAICweights=TRUE,
rmNaParameters=TRUE)

##Save results to a .txt file
write(t(totalData), file = "totalData.txt", ncolumns = 40, append = FALSE, 
sep = "\t")	
}	
}

