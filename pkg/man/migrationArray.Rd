\name{migration.array}
\alias{migrationArray_CollapseForbidden_npop2_maxK1}
\alias{migrationArray_CollapseForbidden_npop2_maxK10}
\alias{migrationArray_CollapseForbidden_npop2_maxK2}
\alias{migrationArray_CollapseForbidden_npop2_maxK3}
\alias{migrationArray_CollapseForbidden_npop2_maxK4}
\alias{migrationArray_CollapseForbidden_npop2_maxK5}
\alias{migrationArray_CollapseForbidden_npop2_maxK6}
\alias{migrationArray_CollapseForbidden_npop2_maxK7}
\alias{migrationArray_CollapseForbidden_npop2_maxK8}
\alias{migrationArray_CollapseForbidden_npop2_maxK9}
\alias{migrationArray_CollapseForbidden_npop3_maxK1}
\alias{migrationArray_CollapseForbidden_npop3_maxK10}
\alias{migrationArray_CollapseForbidden_npop3_maxK2}
\alias{migrationArray_CollapseForbidden_npop3_maxK3}
\alias{migrationArray_CollapseForbidden_npop3_maxK4}
\alias{migrationArray_CollapseForbidden_npop3_maxK5}
\alias{migrationArray_CollapseForbidden_npop3_maxK6}
\alias{migrationArray_CollapseForbidden_npop3_maxK7}
\alias{migrationArray_CollapseForbidden_npop3_maxK8}
\alias{migrationArray_CollapseForbidden_npop3_maxK9}
\alias{migrationArray_CollapseForbidden_npop4_maxK1}
\alias{migrationArray_CollapseForbidden_npop4_maxK10}
\alias{migrationArray_CollapseForbidden_npop4_maxK2}
\alias{migrationArray_CollapseForbidden_npop4_maxK3}
\alias{migrationArray_CollapseForbidden_npop4_maxK4}
\alias{migrationArray_CollapseForbidden_npop4_maxK5}
\alias{migrationArray_CollapseForbidden_npop4_maxK6}
\alias{migrationArray_CollapseForbidden_npop4_maxK7}
\alias{migrationArray_CollapseForbidden_npop4_maxK8}
\alias{migrationArray_CollapseForbidden_npop4_maxK9}
\alias{migrationArray_CollapseForbidden_npop5_maxK1}
\alias{migrationArray_CollapseForbidden_npop5_maxK2}
\alias{migrationArray_CollapseForbidden_npop5_maxK3}
\alias{migrationArray_CollapseForbidden_npop5_maxK4}
\alias{migrationArray_CollapseForbidden_npop5_maxK5}
\alias{migrationArray_CollapseForbidden_npop5_maxK6}
\alias{migrationArray_CollapseForbidden_npop5_maxK7}
\alias{migrationArray_CollapseForbidden_npop5_maxK8}
\alias{migrationArray_CollapseForbidden_npop6_maxK1}
\alias{migrationArray_FullyResolvedCollapse_npop2_maxK10}
\alias{migrationArray_FullyResolvedCollapse_npop2_maxK2}
\alias{migrationArray_FullyResolvedCollapse_npop2_maxK3}
\alias{migrationArray_FullyResolvedCollapse_npop2_maxK4}
\alias{migrationArray_FullyResolvedCollapse_npop2_maxK5}
\alias{migrationArray_FullyResolvedCollapse_npop2_maxK6}
\alias{migrationArray_FullyResolvedCollapse_npop2_maxK7}
\alias{migrationArray_FullyResolvedCollapse_npop2_maxK8}
\alias{migrationArray_FullyResolvedCollapse_npop2_maxK9}
\alias{migrationArray_FullyResolvedCollapse_npop3_maxK10}
\alias{migrationArray_FullyResolvedCollapse_npop3_maxK3}
\alias{migrationArray_FullyResolvedCollapse_npop3_maxK4}
\alias{migrationArray_FullyResolvedCollapse_npop3_maxK5}
\alias{migrationArray_FullyResolvedCollapse_npop3_maxK6}
\alias{migrationArray_FullyResolvedCollapse_npop3_maxK7}
\alias{migrationArray_FullyResolvedCollapse_npop3_maxK8}
\alias{migrationArray_FullyResolvedCollapse_npop3_maxK9}
\alias{migrationArray_FullyResolvedCollapse_npop4_maxK4}
\alias{migrationArray_FullyResolvedCollapse_npop4_maxK5}
\alias{migrationArray_FullyResolvedCollapse_npop4_maxK6}
\alias{migrationArray_npop2_maxK1}
\alias{migrationArray_npop2_maxK10}
\alias{migrationArray_npop2_maxK2}
\alias{migrationArray_npop2_maxK3}
\alias{migrationArray_npop2_maxK4}
\alias{migrationArray_npop2_maxK5}
\alias{migrationArray_npop2_maxK6}
\alias{migrationArray_npop2_maxK7}
\alias{migrationArray_npop2_maxK8}
\alias{migrationArray_npop2_maxK9}
\alias{migrationArray_npop3_maxK1}
\alias{migrationArray_npop3_maxK10}
\alias{migrationArray_npop3_maxK2}
\alias{migrationArray_npop3_maxK3}
\alias{migrationArray_npop3_maxK4}
\alias{migrationArray_npop3_maxK5}
\alias{migrationArray_npop3_maxK6}
\alias{migrationArray_npop3_maxK7}
\alias{migrationArray_npop3_maxK8}
\alias{migrationArray_npop3_maxK9}
\alias{migrationArray_npop4_maxK1}
\alias{migrationArray_npop4_maxK2}
\alias{migrationArray_npop4_maxK3}
\alias{migrationArray_npop4_maxK4}
\alias{migrationArray_npop4_maxK5}
\alias{migrationArray_npop5_maxK1}
\title{Cached migration arrays}
\description{
  This has cached files with pre-computed sets of models. 
  npop is the number of populations at the present. 
  maxK is the maximum number of free parameters.
  CollapseForbidden models are all island models (no tree structure of populations)
  FullyResolvedCollapse models are fully tree models
  Thus, loading migrationArray_FullyResolvedCollapse_npop2_maxK6 means you will have 
  all possible tree structures for your lineages, but that all lineages must collapse
  to one lineage at the root. There are only two lineages, and at most six free parameters.
}
\usage{
data(migrationArray_FullyResolvedCollapse_npop2_maxK6)
}
\format{
  The data are stored as a migrationArray and a maxK.
  Note that a popVector is also included, but you should overwrite it with your own.
}
\examples{
data(migrationArray_FullyResolvedCollapse_npop2_maxK6)
print(length(migrationArray))
}
\keyword{datasets}
