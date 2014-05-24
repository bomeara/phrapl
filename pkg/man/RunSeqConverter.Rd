\name{RunSeqConverter}
\alias{RunSeqConverter}
\title{
RunSeqConverter
}
\description{
This function runs a perl script that takes one or more sequence files in nexus, fasta, or Se-Al formats and converts them to (relax) phylip format. This format is required for RAxML (RunRaxml). 
}
\usage{
RunSeqConverter(seqConvPath = paste(getwd(), "/", sep = ""), inFilePath = paste(getwd(), "/", sep = ""), inputFormat = "nex")
}
\arguments{
  \item{seqConvPath}{
Path to the seqConverter script (default is the current dir).
}
  \item{inFilePath}{
Path to the sequence files (default is current dir).
}
  \item{inputFormat}{
Format of sequence files. Default input format is "nex", but can also be "fasta" or "seal" (Se-Al). The sequence files must end with one of these suffixes (e.g. *.fas, *.nex, *.seal).
}
}
\details{
The script "seqConverter.pl" can be found in "/home_path/phrapl/pkg/extdata/"

This perl script was written by Olaf R.P. Bininda-Emonds. For more information see: http://www.uni-oldenburg.de/ibu/systematik-evolutionsbiologie/programme/


}
\author{
  Brian O'Meara, Bryan Carstens, Nathan Jackson, Ariadna Morales-Garcia

  Maintainer: Brian O'Meara <bomeara@utk.edu>
}
\references{
  O'Meara, B.C., N. Jackson, A. Morales-Garcia, and B. Carstens (2014) Phrapl in prep.
}
\seealso{
\link{RunRaxml}
}
\examples{
## The script "seqConverter.pl" can be found in "/home_path/phrapl/pkg/extdata/"
#inFilePath<-load("/path/sequence_files")

#RunSeqConverter(seqConvPath="/phrapl/pkg/extdata/seqConverter.pl", inFilePath=inFilePath, inputFormat="nex")
}
