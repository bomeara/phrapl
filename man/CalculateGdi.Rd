\name{CalculateGdi}
\alias{CalculateGdi}
\alias{CalculateGDI}
\title{
Calculate the genealogical divergence index (gdi)
}
\description{
This function calculates the genealogical divergence index (gdi) for a pair of taxa based on specified
bi-directional migration rates (M) and a divergence time (t) estimated from a PHRAPL analysis.
This index provides something similar to the genealogical sorting index (gsi), except for a given set of
migration and divergence values instead of for a group on a tree. It also differs from the gsi in that
it is calculated for a given pair of populations or taxa rather than for a single taxon in respect to the
entire tree. The index is scaled to be between 0 and 1, where 0 = panmixia and 1 = strong divergence.
Except when gdi is close to 1, there is a fairly high error, such that precision requires a large number
of replicates (e.g., nreps = 10,000). We have incorporated the \code{binom.confint} function from the
\code{binom package} to calculate a confidence interval if desired.

}
\usage{
CalculateGdi(tau,migration.in,migration.out=NULL,nreps=10000,ciMethod="exact", msPath =
                    system.file("msdir", "ms", package = "P2C2M"))
}
\arguments{
  \item{tau}{
The divergence time between two populations/groups.
}
  \item{migration.in}{
The migration rate into a population
}
  \item{migration.out}{
The migration rate out of a population. If this is not specified, migration rates are assumed
to be symmetrical.
}
  \item{nreps}{
The number of genealogies to be simulated.
}
  \item{ciMethod}{
Which method to use to calculate a confidence interval for gdi. Any or all of eight methods
for calculating a binomial confidence interval can be specified ("exact", "asymptotic",
"agresti-coull", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit", or "profile").
\code{ciMethod = "all"} will result in all eight being calculated. The default method is the
Clopper-Pearson interval ("exact"). If \code{ciMethod=NULL}, only the gdi will be outputted.
}
\item{msPath}{
Path to the local installation of ms; typing this string on the command line
should result in ms running.
}
}
\author{
  Nathan Jackson and Brian O'Meara
}
\examples{
#CalculateGdi(tau=1,migration.in=1)
}
