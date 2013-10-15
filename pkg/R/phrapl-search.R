#has the functions for searching or optimizing a single model

# searchContinuousModelSpace<-function(p, migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, maxParameterValue=100, nTrees=1 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=FALSE, debug=FALSE,method="nlminb",itnmax=NULL,...) {
  # modelID<-returnModel(p,migrationArrayMap)
  # if(print.results) {
    # resultVector<-c(modelID,p)
    # names(resultVector)<-c("migrationArryIndividualID","collapseMatrix.number", "n0multiplierMap.number","migrationArray.number")
    # print(resultVector)
  # }
  # if(is.na(modelID)) {
    # return(badAIC)
  # }
  # else {
    # paramNames<-msIndividualParameters(migrationArray[[modelID]])
    # startingVals<-log(c(rlnorm(sum(grepl("collapse",paramNames)),1,1), rlnorm(sum(grepl("n0multiplier",paramNames)),1,1), rbeta(sum(grepl("migration",paramNames)),shape1=1,shape2=3) )) #going to optimize in log space
    # if(debug) {
      # print(startingVals) 
    # }
    # searchResults<-optimx(par=startingVals, fn=eeturnAIC, method=method, migrationIndividual=migrationArray[[modelID]], migrationArrayMap=migrationArrayMap, migrationArray=migrationArray, popVector=popVector, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, print.ms.string=print.ms.string, print.results=print.results, itnmax=itnmax, debug=debug, control=list(iter.max=itnmax, maxit=itnmax,iter.lim=itnmax),...)
    # return(searchResults$fvalues)
  # }
# }

SearchContinuousModelSpaceOptim<-function(p, migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, maxParameterValue=100, nTrees=1 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=FALSE, debug=FALSE,method="nlminb",itnmax=NULL, return.all=FALSE, ...) {
  modelID<-ReturnModel(p,migrationArrayMap)
  if(print.results) {
    resultVector<-c(modelID,p)
    names(resultVector)<-c("migrationArryIndividualID","collapseMatrix.number", "n0multiplierMap.number","migrationArray.number")
    print(resultVector)
  }
  if(is.na(modelID)) {
    return(badAIC)
  }
  else {
    paramNames<-MsIndividualParameters(migrationArray[[modelID]])
    startingVals<-log(c(rlnorm(sum(grepl("collapse",paramNames)),1,1), rlnorm(sum(grepl("n0multiplier",paramNames)),1,1), rbeta(sum(grepl("migration",paramNames)),shape1=1,shape2=3) )) #going to optimize in log space
    if(debug) {
      print(startingVals) 
    }
    searchResults<-optim(par=startingVals, fn=ReturnAIC, method=method, control=list(maxit=itnmax), migrationIndividual=migrationArray[[modelID]], popVector=popVector, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, print.ms.string=print.ms.string, print.results=print.results, debug=debug, ...)
    ifelse(return.all, return(searchResults), return(searchResults$value))     
  }
}

SearchContinuousModelSpaceNLoptr<-function(p, migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, maxParameterValue=100, nTrees=1 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=FALSE, debug=FALSE,method="nlminb",itnmax=NULL, return.all=FALSE, maxtime=0, maxeval=0, ...) {
  modelID<-ReturnModel(p,migrationArrayMap)
  if(print.results) {
    resultVector<-c(modelID,p)
    names(resultVector)<-c("migrationArryIndividualID","collapseMatrix.number", "n0multiplierMap.number","migrationArray.number")
    print(resultVector)
  }
  if(is.na(modelID)) {
    return(badAIC)
  }
  else {
    paramNames<-MsIndividualParameters(migrationArray[[modelID]])
    startingVals<-log(c(rlnorm(sum(grepl("collapse",paramNames)),1,1), rlnorm(sum(grepl("n0multiplier",paramNames)),1,1), rbeta(sum(grepl("migration",paramNames)),shape1=1,shape2=3) )) #going to optimize in log space
    if(debug) {
      print(startingVals) 
    }
    searchResults<-nloptr(x0=startingVals, eval_f=ReturnAIC, opts=list("maxeval"=itnmax, "algorithm"="NLOPT_LN_SBPLX", "print_level"=1, maxtime=maxtime, maxeval=maxeval), migrationIndividual=migrationArray[[modelID]], popVector=popVector, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, print.ms.string=print.ms.string, print.results=print.results, debug=debug)
    ifelse(return.all, return(searchResults), return(searchResults$objective))     
  }
}


# searchDiscreteModelSpace<-function(migrationArrayMap, migrationArray, popVector, print.ms.string=FALSE, badAIC=100000000000000, maxParameterValue=100, nTrees=1,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, debug=FALSE,  method="newuoa", itnmax=NULL,pop.size=50, print.results=FALSE, ...) {
  # Domains<-matrix(ncol=2,nrow=3)
  # Domains[1,]<-range(migrationArrayMap$collapseMatrix.number)
  # Domains[2,]<-range(migrationArrayMap$n0multiplierMap.number)
  # Domains[3,]<-range(migrationArrayMap$migrationArray.number)
  
  # results<-genoud(searchContinuousModelSpace,nvars=3, max=FALSE,starting.values=c(1,1,1), MemoryMatrix=TRUE, boundary.enforcement=2, data.type.int=TRUE, Domains=Domains, migrationArrayMap=migrationArrayMap, migrationArray=migrationArray, popVector=popVector, print.ms.string=print.ms.string, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, debug=debug, method=method,itnmax=itnmax, pop.size=pop.size, print.results=print.results,...)
  # return(results)
# }

SearchDiscreteModelSpaceOptim<-function(migrationArrayMap, migrationArray, popVector, print.ms.string=FALSE, badAIC=100000000000000, maxParameterValue=100, nTrees=1,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, debug=FALSE,  method="BFGS", itnmax=NULL,pop.size=50, print.results=FALSE, ...) {
  Domains<-matrix(ncol=2,nrow=3)
  Domains[1,]<-range(migrationArrayMap$collapseMatrix.number)
  Domains[2,]<-range(migrationArrayMap$n0multiplierMap.number)
  Domains[3,]<-range(migrationArrayMap$migrationArray.number)
  
  results<-genoud(SearchContinuousModelSpaceOptim,nvars=3, max=FALSE,starting.values=c(1,1,1), MemoryMatrix=TRUE, boundary.enforcement=2, data.type.int=TRUE, Domains=Domains, migrationArrayMap=migrationArrayMap, migrationArray=migrationArray, popVector=popVector, print.ms.string=print.ms.string, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, debug=debug, method=method,itnmax=itnmax, pop.size=pop.size, print.results=print.results,...)
  return(results)
}

SearchDiscreteModelSpaceNLoptr<-function(migrationArrayMap, migrationArray, popVector, print.ms.string=FALSE, badAIC=100000000000000, maxParameterValue=100, nTrees=1,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, debug=FALSE,  method="nlminb", itnmax=NULL,pop.size=50, print.results=FALSE, maxtime=0, maxeval=0, ...) {
  Domains<-matrix(ncol=2,nrow=3)
  Domains[1,]<-range(migrationArrayMap$collapseMatrix.number)
  Domains[2,]<-range(migrationArrayMap$n0multiplierMap.number)
  Domains[3,]<-range(migrationArrayMap$migrationArray.number)
  
  results<-genoud(SearchContinuousModelSpaceNLoptr,nvars=3, max=FALSE,starting.values=c(1,1,1), MemoryMatrix=TRUE, boundary.enforcement=2, data.type.int=TRUE, Domains=Domains, migrationArrayMap=migrationArrayMap, migrationArray=migrationArray, popVector=popVector, print.ms.string=print.ms.string, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees,msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, debug=debug, method=method,itnmax=itnmax, pop.size=pop.size, print.results=print.results, maxtime=maxtime, maxeval=maxeval, return.all=FALSE,...)
  return(results)
}



ExhaustiveSearchOptim<-function(migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, maxParameterValue=100, nTrees=1 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=FALSE, debug=FALSE,method="nlminb",itnmax=NULL, ncores=1, results.file=NULL, return.all=TRUE, ...) {
  AIC.values<-rep(NA,length(migrationArray))
  for (i in sequence(length(migrationArray))) {
  	p<-c(migrationArrayMap$collapseMatrix.number[i], migrationArrayMap$n0multiplierMap.number[i], migrationArrayMap$migrationArray.number[i])
  	try(AIC.values[i]<-SearchContinuousModelSpaceOptim(p, migrationArrayMap, migrationArray, popVector, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees, msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, print.ms.string=print.ms.string, print.results=print.results, debug=debug,method=method,itnmax=itnmax,...))
  	print(c(i, length(migrationArray), i/length(migrationArray), AIC.values[i]))
	if(!is.null(results.file)) {
		save(list=ls(), file=results.file)
	}
  }
  return(AIC.values)
}

ExhaustiveSearchNLoptr<-function(migrationArrayMap, migrationArray, popVector, badAIC=100000000000000, maxParameterValue=100, nTrees=1 ,msLocation="/usr/local/bin/ms",compareLocation="comparecladespipe.pl",assign="assign.txt",observed="observed.txt",unresolvedTest=TRUE, print.ms.string=FALSE, print.results=FALSE, debug=FALSE,method="nlminb",itnmax=NULL, ncores=1, results.file=NULL, maxtime=0, maxeval=0, return.all=TRUE, ...) {
  AIC.values<-rep(NA,length(migrationArray))
  results.list<-list()
  for (i in sequence(length(migrationArray))) {
  	p<-c(migrationArrayMap$collapseMatrix.number[i], migrationArrayMap$n0multiplierMap.number[i], migrationArrayMap$migrationArray.number[i])
  	if(return.all) {
  		result.indiv<-NULL
  		try(result.indiv<-SearchContinuousModelSpaceNLoptr(p, migrationArrayMap, migrationArray, popVector, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees, msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, print.ms.string=print.ms.string, print.results=print.results, debug=debug,method=method,itnmax=itnmax, maxtime=maxtime, maxeval=maxeval, return.all=return.all, ...))
  		print(result.indiv)
  		if(!is.null(result.indiv)) {
  			AIC.values[i]<-result.indiv$objective	
  		}
  		results.list<-append(results.list, list(result.indiv))
  	} else {
  		try(AIC.values[i]<-SearchContinuousModelSpaceNLoptr(p, migrationArrayMap, migrationArray, popVector, badAIC=badAIC, maxParameterValue=maxParameterValue, nTrees=nTrees, msLocation=msLocation,compareLocation=compareLocation,assign=assign,observed=observed,unresolvedTest=unresolvedTest, print.ms.string=print.ms.string, print.results=print.results, debug=debug,method=method,itnmax=itnmax, maxtime=maxtime, maxeval=maxeval, return.all=return.all, ...))
  	}
  	 print(c(i, length(migrationArray), i/length(migrationArray), AIC.values[i]))

  	if(!is.null(results.file)) {
		save(list=ls(), file=results.file)
	}

  }
  ifelse(return.all, return(results.list), return(AIC.values))     
}
