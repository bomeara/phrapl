\name{CompileMS}
\alias{CompileMS}
\alias{CompileMs}
\title{
Compile ms
}
\description{
PHRAPL uses the program ms (Hudson 2002) to simulate gene trees and we have packaged the raw C code for the 
program within PHRAPL. Thus, upon downloading PHRAPL, users should compile ms using the \code{CompileMS}
function. To check whether the compilation was successful, you can run the function \code{CheckMS}.

}
\usage{
CheckMS(rand=1)
}
\arguments{
  \item{rand}{
ms can use one of two different random number generators. If compilation using the default \code{rand = 1} 
does not work, try \code{rand = 2}.  
}
}
\references{
Hudson, R.R. 2002. Generating samples under a Wright-Fisher neutral model. \emph{Bioinformatics}. 18:337-338.
\url{http://home.uchicago.edu/~rhudson1/source.html}
}
\author{
  Nathan Jackson
}
