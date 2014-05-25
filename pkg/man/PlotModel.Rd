\name{PlotModel}
\alias{PlotModel}
\alias{plotmodel}
\alias{Plotmodel}
\alias{plot.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a 3D spinning plot for a given model
}
\description{
This uses rgl to create a depiction of a given model, complete with migration 
arrows. You can spin this using the cursor.
}
\usage{
PlotModel(migrationIndividual, parameterVector = NULL, taxonNames = NULL, 
time.axis = FALSE, time.axis.col = "black", apply.base = FALSE, base.color = 
"black", new.window = TRUE)
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
  \item{time.axis}{
If you want to have a time axis, enter TRUE.
}
  \item{time.axis.col}{
If you want a particular color for the time axis, enter it here.
}
  \item{apply.base}{
This will draw a box at the rootward end of the population history. This can be 
convenient for 3D printing.
}
  \item{base.color}{
What color should the box be, if one is drawn.
}
  \item{new.window}{
If FALSE, overwrites the current image if one exists. By default, TRUE, so each 
time you run this function you create another window.
}
}
\details{
This requires rgl, which periodically becomes hard to install on certain 
platforms.
}
\value{
Creates a stunning image.
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
