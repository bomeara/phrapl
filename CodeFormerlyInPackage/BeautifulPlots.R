library(phrapl)
setwd("~/Dropbox/PhraplExample/PostRunAnalysis")
data<-read.delim("~/Dropbox/PhraplExample/PostRunAnalysis/results.all_2014-05-19-empirical.txt", stringsAsFactors=FALSE)
load("~/Dropbox/PhraplExample/PostRunAnalysis/MigrationArray_3pop_4K.rda")
migration.model <- 4
speciestree.models <- c(65, 21, 43)
datasets <- unique(data$dataset)
datasets.bestwt<-rep(NA, length(datasets))
for(dataset.index in sequence(length(datasets))) {
		data.subset<-data[which(data$dataset==datasets$Data[dataset.index]),]
		datasets.bestwt[dataset.index]<-max(data.subset$wAIC)
}
datasets<-datasets[order(datasets.bestwt)]


GetDistanceFromParameterEstimates<-function(data.subset, models, na.replace=0) {
	rows<-match(models, data.subset$models)
	parameters<-data.subset[rows, (grepl("_", colnames(data), perl=FALSE) * !grepl("_I", colnames(data), perl=FALSE))]	
	parameters[is.na(parameters)]<-na.replace
	return(as.numeric(dist(parameters)))
}

GetDistanceFromParameterIdentities<-function(data.subset, models, na.replace=0.5) {
	rows<-match(models, data.subset$models)
	parameters<-data.subset[rows, grepl("_I", colnames(data), perl=FALSE)]	
	parameters[is.na(parameters)]<-na.replace
	return(as.numeric(dist(parameters)))
}

GetMigrationDistance<-function(focal.model, data.subset, migration.model=migration.model, criterion="estimates") {
	if(criterion=="estimates") {
		return(GetDistanceFromParameterEstimates(data.subset, c(focal.model, migration.model)))
	}		
	if(criterion=="identities") {
		return(GetDistanceFromParameterIdentities(data.subset, c(focal.model, migration.model)))
	}		
}

GetSpeciesDistance<-function(focal.model, data.subset, speciestree.models= speciestree.models, criterion="estimates") {
	if(criterion=="estimates") {
		GetDistance <- GetDistanceFromParameterEstimates
	}		
	if(criterion=="identities") {
		GetDistance <- 	GetDistanceFromParameterIdentities
	}
	min.dist<-Inf
	for (i in sequence(length(speciestree.models))) {
		min.dist<-min(min.dist, GetDistance(data.subset, models=c(focal.model, speciestree.models[i])))	
	}		
	return(min.dist)
}




#data.subset<-data[which(data$dataset==datasets[1]),]
#data.subset.distances <- data.frame(data.subset, migrationDistanceParams = sapply(unique(data$models), GetMigrationDistance, data.subset=data.subset, migration.model = migration.model), speciesDistanceParams = sapply(unique(data$models), GetSpeciesDistance, data.subset=data.subset, speciestree.models = speciestree.models), migrationDistanceIdent = sapply(unique(data$models), GetMigrationDistance, data.subset=data.subset, migration.model = migration.model, criterion="identities"), speciesDistanceIdent = sapply(unique(data$models), GetSpeciesDistance, data.subset=data.subset, speciestree.models = speciestree.models, criterion="identities"))
#data.subset.distances<-data.subset.distances[order(data.subset$models),]
#model.order<-order(data.subset.distances$migrationDistanceIdent, decreasing=TRUE)

KCollapse <- function(migrationIndividual) {
	return(KCollapseMatrix(migrationIndividual$collapseMatrix))	
}



models.fractionNonZeroMigration <- sapply(migrationArray,FractionNonZeroMigration)
models.KCollapse <- sapply(migrationArray, KCollapse)
models.NumberPopulationsAtRoot <- sapply(migrationArray, NumberPopulationsAtRoot)
ComputeModelOrder <- function(fnzm, kc, np) { #want decreasing=FALSE, so first ones are species models and have lowest values, last ones are migration only
	if(kc==0) { #migration only
		return(50+fnzm)
	}	
	if(fnzm==0) { #collapse only
		return(-.01*kc)
	}
	if(kc==2) { #tree
		return(2 + fnzm)
	}
	if (kc==1) { #could be tree, or just collapse of one
		if(np==1) { #star tree
			return(4 + fnzm)
		} else {
			return(6 + fnzm)
		}
	}
	stop("all should be handled")
}
models.scores <- rep(NA, length(models.KCollapse)) 

for (i in sequence(length(models.scores))) {
	models.scores[i] <- ComputeModelOrder(models.fractionNonZeroMigration[i], models.KCollapse[i], models.NumberPopulationsAtRoot[i])
}
model.order<-order(models.scores)

#model.order<-order(sapply(migrationArray,FractionNonZeroMigration), -sapply(migrationArray, KCollapse))

plot(x=c(-11, length(model.order)+1), y=c(0, length(datasets)+5), bty="n", xlab="", xaxt="n", type="n", yaxt="n", ylab="")
polygon(x=c(0.5,0.5,4.5,4.5), y=c(.5,.5+length(datasets),.5+length(datasets),.5), col="lightgray", border=NA)
polygon(x=c(length(model.order)+0.5,length(model.order)+0.5,length(model.order)-3.5,length(model.order)-3.5), y=c(.5,.5+length(datasets),.5+length(datasets),.5), col="lightgray", border=NA)


text(x=0, y=sequence(length(datasets)), labels=gsub("_", " & ",datasets), pos=2, cex=0.5)
text(x=2.5, y=0, labels="Tree")
text(x=length(model.order)-1.5, y=0, labels="Island")
#axis(side=1, at=c(1, length(model.order)), labels=c("Tree, no migration", "Island model"))
for(dataset.index in sequence(length(datasets))) {
	data.subset<-data[which(data$dataset==datasets[dataset.index]),]
	data.subset.distances <- data.frame(data.subset, migrationDistanceParams = sapply(unique(data.subset$models), GetMigrationDistance, data.subset=data.subset, migration.model = migration.model), speciesDistanceParams = sapply(unique(data.subset$models), GetSpeciesDistance, data.subset=data.subset, speciestree.models = speciestree.models), migrationDistanceIdent = sapply(unique(data.subset$models), GetMigrationDistance, data.subset=data.subset, migration.model = migration.model, criterion="identities"), speciesDistanceIdent = sapply(unique(data.subset$models), GetSpeciesDistance, data.subset=data.subset, speciestree.models = speciestree.models, criterion="identities"))
data.subset.distances<-data.subset.distances[order(data.subset$models),]
	lines(x=c(1, length(model.order)), y=rep(dataset.index, 2), col="gray", lwd=0.5)
	for (i in sequence(dim(data.subset.distances)[1])) {
		matching.model <- which(data.subset.distances$models==model.order[i])
		if(length(matching.model)==1) {
			col="black"
			if(data.subset.distances$wAIC[matching.model]==max(data.subset.distances$wAIC) && data.subset.distances$wAIC[matching.model] > 0.5) {
				col="red"
			}
			symbols(i, dataset.index, circles=sqrt(data.subset.distances$wAIC[matching.model])*sqrt(1/pi), bg=col, fg=NULL, add=TRUE, inches=FALSE)	#circle area scales with the weight, not diameter. 
		} else {
			print(paste("no match for ", model.order[i], "in", data.subset.distances$dataset[1]))	
		}
	}
	
}
#text(0, length(datasets)+2, "K collapse", cex=0.7, pos=2)
#text(0, length(datasets)+3, "Fraction migration", cex=0.7, pos=2)
#text(0, length(datasets)+4, "Number pops at root", cex=0.7, pos=2)

for (i in sequence(length(migrationArray))) {
#	text(i, length(datasets)+2, KCollapseMatrix(migrationArray[[model.order[i]]]$collapseMatrix), cex=0.5)	
}
for (i in sequence(length(migrationArray))) {
#	text(i, length(datasets)+3, round(FractionNonZeroMigration(migrationArray[[model.order[i]]]),1), cex=0.5)	
}
for (i in sequence(length(migrationArray))) {
#	text(i, length(datasets)+4, round(NumberPopulationsAtRoot(migrationArray[[model.order[i]]]),1), cex=0.5)	
}

