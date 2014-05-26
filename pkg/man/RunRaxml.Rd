\name{RunRaxml}
\alias{RunRaxml}
\title{
RunRaxml
}
\description{
This function runs RAxML through phrapl. The function will produce an input 
string for RAxML and calling it up. RAxML must be installed in the system before
invoking this function. 
It reads in all phylip files in the designated path, and runs RAxML for each in
turn. Outgroups and mutation models can be specified either as a single string
to be used for all loci or as a vector which needs to match the order of the 
reading of the phylip files (i.e., alphabetic/numeric).
}
\usage{
RunRaxml(raxmlPath = "./", raxmlVersion = "raxmlHPC",inputPath = "./", mutationModel, 
outgroup, iterations, seed = sample(1:1e+07, 1), outputSeeds = FALSE, discard = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raxmlPath}{
Path to the RAxML executable (default is current dir).
}
  \item{raxmlVersion}{
Name and version of the executable (default is "raxmlHPC").
}
  \item{inputPath}{
Path to sequence files in phylip format.
}
  \item{mutationModel}{
Mutation model (e.g., "GTRGAMMA"). If there is only one, it can be simply given
as a single string. 
If there is a different model for each locus, the path to a text file can be
provided that includes a single column of models that is the same length as the
number of sequence files (note that the ordering must match the order of the 
phylip files in the directory (i.e., alphabetic/numeric)).
}
  \item{outgroup}{
A character string denoting an outgroup. Like with mutation models, a single 
outgroup can be given to use for all sequence files, or the location of a text 
file can be given which provides a vector (a single column) with the outgroups.
}
  \item{iterations}{
The number of iterations of the search algorithm.
}
  \item{seed}{
Provide a seed (a number). The default is a random number
}
  \item{outputSeeds}{
If "outputSeeds=TRUE", a file is saved that includes the seeds used for each 
dataset
}
  \item{discard}{
If "discard=TRUE", all RAxML output except the most likely tree is automatically 
discarded
}
}
\details{
For more information see http://bodegaphylo.wikispot.org/RAxML_Tutorial.
}
\author{
  Brian O'Meara, Bryan Carstens, Nathan Jackson, Ariadna Morales-Garcia

  Maintainer: Brian O'Meara <bomeara@utk.edu>
}
\references{
  O'Meara, B.C., N. Jackson, A. Morales-Garcia, and B. Carstens (2014) 
  Phrapl in prep.
}
\seealso{
\link{RunSeqConverter}
}
\examples{
#raxmlPath<-load("/path_to_RAxML_executable/")
#inputPath<-load("/path/sequence_files.phylip")

#RunRaxml(raxmlPath=raxmlPath, 
#raxmlVersion="raxmlHPC",
#inputPath=inputPath,
#mutationModel="GTRGAMMA",
#outgroup="outgroup",
#iterations=2,
#seed=sample(1:10000000,1),
#outputSeeds=FALSE,
#discard=FALSE)

## The output file will be located in your current dir and will be named 
##"RAxML_bestTree.YOUR_FILE_NAME.phylip.tre"	

}
