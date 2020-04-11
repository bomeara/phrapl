# test_that("GridSearch works with P2C2M install", {
#   data(TestData)
# observedTrees<-PrepSubsampling(assignmentsGlobal=assignFile,
#   observedTrees<-trees,popAssignments=popAssignments,
#   subsamplesPerGene=10)
# subsampleWeights.df<-GetPermutationWeightsAcrossSubsamples(
#   popAssignments=popAssignments,observedTrees=observedTrees)
# migrationArrayMap <- GenerateMigrationArrayMap(migrationArray)
# result<-GridSearch(migrationArrayMap=migrationArrayMap,migrationArray=
#   migrationArray,modelRange=1:6,popAssignments=popAssignments,
#   observedTrees=observedTrees,subsampleWeights.df=subsampleWeights.df,
#   subsamplesPerGene=10,totalPopVector=rep(20,10),nTrees=1e3,
#   print.results=TRUE,return.all=TRUE, usePhyclust=FALSE)
#   expect_lte(max(result$overall.results$lnL),-50)
# })

test_that("GridSearch works with phyclust install", {
  data(TestData)
observedTrees<-PrepSubsampling(assignmentsGlobal=assignFile,
  observedTrees<-trees,popAssignments=popAssignments,
  subsamplesPerGene=10)
subsampleWeights.df<-GetPermutationWeightsAcrossSubsamples(
  popAssignments=popAssignments,observedTrees=observedTrees)
migrationArrayMap <- GenerateMigrationArrayMap(migrationArray)
result<-GridSearch(migrationArrayMap=migrationArrayMap,migrationArray=
  migrationArray,modelRange=1:6,popAssignments=popAssignments,
  observedTrees=observedTrees,subsampleWeights.df=subsampleWeights.df,
  subsamplesPerGene=10,totalPopVector=rep(20,10),nTrees=1e3,
  print.results=TRUE,return.all=TRUE, usePhyclust=TRUE)
  expect_lte(max(result$overall.results$lnL),-50)
})
