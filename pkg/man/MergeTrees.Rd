\name{MergeTrees}
\alias{MergeTrees}
\title{
MergeTrees
}
\description{
This function merges all the *.tre files in a directory into a single file called "trees.tre". The files must be in newick format. These tree files must have the same label for individuals across all loci. These labels must be the same listed in assignFile ("cladeAssignments.txt").
}
\usage{
MergeTrees(treesPath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treesPath}{
Path to tree files in newick format.
}
}
\author{
  Brian O'Meara, Bryan Carstens, Nathan Jackson, Ariadna Morales-Garcia

  Maintainer: Brian O'Meara <bomeara@utk.edu>
}
\references{
  O'Meara, B.C., N. Jackson, A. Morales-Garcia, and B. Carstens (2014) Phrapl in prep.
}

\examples{

#MergeTrees(treePath="/path_to_tree_files/")

}
