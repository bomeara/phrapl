\name{phrapl-data}
\alias{CachedMigrationArrays}
\alias{MigrationArray_3pops_maxMigrationK1_asymmetrical}
\alias{MigrationArray_3pops_maxMigrationK1}
\alias{MigrationArray_3pops_maxMigrationK2_asymmetrical}
\alias{MigrationArray_3pops_maxMigrationK2}
\alias{MigrationArray_3pops_maxMigrationK3_asymmetrical}
\alias{MigrationArray_3pops_maxMigrationK3}
\alias{MigrationArray_3pops_maxMigrationK4_asymmetrical}
\alias{MigrationArray_4pops_maxMigrationK1}
\alias{MigrationArray_4pops_maxMigrationK2}

\title{
Cached migrationArrays
}
\description{
These are precomputed migrationArrays that have been generated using
GenerateMigrationIndividuals according to specific critera. All the model
sets avaiable below contain only models with two type of parameters: 
coalescence times and migration rates (i.e., equal and stable population 
sizes are invoked in all models) and fully resolved trees (i.e., forceTree = TRUE).

To load a migrationArray with 3 populations, a total maximum of three parameters 
(maxK = 3), a maximum of one migration rate (MaxMigrationK = 1), and the 
possibility of assymetrical migration (forceSymmetricalMigration = FALSE), type

\code{data(MigrationArray_3pops_maxMigrationK1_asymmetrical)}.

To load a migrationArray with three populations, a total maximum of three parameters 
(maxK = 3), a maximum of one migration rate (MaxMigrationK = 1), and only
symmetrical migration (forceSymmetricalMigration = TRUE), type

\code{data(MigrationArray_3pops_maxMigrationK1)}.

To load a migrationArray with three populations, a total maximum of four parameters 
(maxK = 4), a maximum of two migration rates (MaxMigrationK = 2), and the 
possibility of assymetrical migration (forceSymmetricalMigration = FALSE), type

\code{data(MigrationArray_3pops_maxMigrationK2_asymmetrical)}.

To load a migrationArray with three populations, a total maximum of four parameters 
(maxK = 4), a maximum of two migration rates (MaxMigrationK = 2), and only
symmetrical migration (forceSymmetricalMigration = TRUE), type

\code{data(MigrationArray_3pops_maxMigrationK2)}.

To load a migrationArray with three populations, a total maximum of five parameters 
(maxK = 5), a maximum of three migration rates (MaxMigrationK = 3), and the 
possibility of assymetrical migration (forceSymmetricalMigration = FALSE), type

\code{data(MigrationArray_3pops_maxMigrationK3_asymmetrical)}.

To load a migrationArray with three populations, a total maximum of five parameters 
(maxK = 5), a maximum of three migration rates (MaxMigrationK = 3), and only
symmetrical migration (forceSymmetricalMigration = TRUE), type

\code{data(MigrationArray_3pops_maxMigrationK3)}.

To load a migrationArray with three populations, a total maximum of six parameters 
(maxK = 6), a maximum of four migration rates (MaxMigrationK = 4), and the 
possibility of assymetrical migration (forceSymmetricalMigration = FALSE), type

\code{data(MigrationArray_3pops_maxMigrationK4_asymmetrical)}.

To load a migrationArray with four populations, a total maximum of four parameters 
(maxK = 4), a maximum of one migration rate (MaxMigrationK = 1), and only
symmetrical migration (forceSymmetricalMigration = TRUE), type

\code{data(MigrationArray_4pops_maxMigrationK1)}.

To load a migrationArray with four populations, a total maximum of five parameters 
(maxK = 5), a maximum of two migration rates (MaxMigrationK = 2), and only
symmetrical migration (forceSymmetricalMigration = TRUE), type

\code{data(MigrationArray_4pops_maxMigrationK2)}.


}
\examples{
  #Load migrationArray		
  data(MigrationArray_3pops_maxMigrationK1_asymmetrical)
  #Inspect the number of models
  length(migrationArray)
}
\keyword{example}
\keyword{TestData}
