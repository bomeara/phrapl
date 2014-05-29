\name{SaveMovie}
\alias{SaveMovie}
\alias{savemovie}
\alias{Savemovie}
\alias{save.movie}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function converts a 3D tree plot into an animated gif, suitable for 
inserting into a slide presentation, putting on a website, etc.
}
\description{
Uses the animation package to save a series of images as a population history 
rotates, then stitches them to form an animated gif.
}
\usage{
SaveMovie(total.revolutions = 1, duration = 10, save.dir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total.revolutions}{
How many full rotations to take.
}
  \item{duration}{
How long should the movie be, in seconds?
}
  \item{save.dir}{
Where should the final movie be saved?
}
}
\details{
This will require the animation package as well as ImageMagick or 
GraphicsMagick. It could be modified to create Quicktime or Flash movies 
instead. 
}
\value{
This will create an animated gif in the specified folder.
}
\author{
Brian O'Meara
}
\seealso{
\code{\link{PlotModel}}, ~~~
}
\keyword{ ~plot }
\keyword{ ~movie }