\name{PlotModel2D}
\alias{PlotModel2D}
\alias{plotmodel2d}
\alias{Plotmodel2d}
\alias{plot.model2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a 2D plot for a given model
}
\description{
This creates a two-dimensional model
}
\usage{
PlotModel2D(migrationIndividual, parameterVector=NULL, taxonNames=NULL, 
	pad.for.input=FALSE, tree.lwd=30, tree.col="black", arrow.lwd=1, 
	arrow.col="red", arrow.lty="solid", tip.cex = 1, tip.col = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{migrationIndividual}{
This is the chosen model. If you want model 17 from migrationArray, you could 
pass this in as migrationIndividual=migrationArray[[17]]
}
  \item{parameterVector}{
If you have numeric parameter estimates, you can pass them in here. This affects
things like node height, migration rate arrow width, etc. Otherwise, it will 
plot a generic model.
}
  \item{taxonNames}{
If instead of population A, B, ... you want Oahu, Hawaii, ... you can pass them 
in as a vector.
}
  \item{pad.for.input}{
[not implemented yet]
}
  \item{tree.lwd}{
Line width for the tree
}
  \item{tree.col}{
Color for the tree
}
  \item{arrow.lwd}{
Line width for migration arrows
}
  \item{arrow.col}{
Color for migration arrows
}
  \item{arrow.lty}{
Line type for migration arrows
}
  \item{tip.cex}{
Specify tip size
}
  \item{tip.col}{
Specify tip color
}
}
\details{
This will create a plot showing migration, population collapse, and other features. It does not yet use parameter estimates. It shows a gradient for branches at the base of the tree to represent the history being lost in the mists of time.
}
\value{
Creates a stunning image.
}
\examples{
data(MigrationArray_3pops_maxMigrationK1)
pretty.model <- migrationArray[[1]]
PlotModel2D(pretty.model)	
}
\author{
Brian O'Meara
}

\seealso{
\code{\link{SaveMovie}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~3D }
\keyword{ ~plot }
