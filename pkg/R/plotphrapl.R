library(rgl)
library(phrapl)

load("/Users/bomeara/Documents/MyDocuments/Active/phrapl/pkg/data/migrationArray_npop3_maxK5.Rsave")


arrow3d2 <- function(base, tip, rad=1, col="red") {
	new.pos<-cbind(seq(from=base[1], to=tip[1], length.out=10), seq(from=base[2], to=tip[2], length.out=10), seq(from=base[3], to=tip[3], length.out=10))
                   	#base<-new.pos[6,]
                   	#tip<-new.pos[8,]
                   	lines3d(x=new.pos[1:8,1], y=new.pos[1:8,2], z=new.pos[1:8,3], col=col, lwd=rad*20)

		shade3d(addNormals(subdivision3d(cylinder3d(rbind(new.pos[8,], new.pos[6,]), radius=c(0.00001, rad), sides=8, twist=1, closed=-2), depth=2)), col=col)
}


#based on cone3d from RGL demo
arrow3d <- function(base=c(0,0,0),tip=c(0,0,1),rad=1,n=30,draw.base=TRUE,qmesh=FALSE,
                   trans = par3d("userMatrix"), ...) {
                   	bounds<-rbind(base,tip)
                   	new.pos<-cbind(seq(from=base[1], to=tip[1], length.out=10), seq(from=base[2], to=tip[2], length.out=10), seq(from=base[3], to=tip[3], length.out=10))
                   	base<-new.pos[6,]
                   	tip<-new.pos[8,]
                   	lines3d(x=new.pos[1:8,1], y=new.pos[1:8,2], z=new.pos[1:8,3])
  ax <- tip-base
  if (missing(trans) && !rgl.cur()) trans <- diag(4)
  ### is there a better way?
  if (ax[1]!=0) {
    p1 <- c(-ax[2]/ax[1],1,0)
    p1 <- p1/sqrt(sum(p1^2))
    if (p1[1]!=0) {
      p2 <- c(-p1[2]/p1[1],1,0)
      p2[3] <- -sum(p2*ax)
      p2 <- p2/sqrt(sum(p2^2))
    } else {
      p2 <- c(0,0,1)
    }
  } else if (ax[2]!=0) {
    p1 <- c(0,-ax[3]/ax[2],1)
    p1 <- p1/sqrt(sum(p1^2))
    if (p1[1]!=0) {
      p2 <- c(0,-p1[3]/p1[2],1)
      p2[3] <- -sum(p2*ax)
      p2 <- p2/sqrt(sum(p2^2))
    } else {
      p2 <- c(1,0,0)
    }
  } else {
    p1 <- c(0,1,0); p2 <- c(1,0,0)
  }
  degvec <- seq(0,2*pi,length=n+1)[-1]
  ecoord2 <- function(theta) {
    base+rad*(cos(theta)*p1+sin(theta)*p2)
  }
  i <- rbind(1:n,c(2:n,1),rep(n+1,n))
  v <- cbind(sapply(degvec,ecoord2),tip)
  if (qmesh) 
    ## minor kluge for quads -- draw tip twice
    i <- rbind(i,rep(n+1,n))
  if (draw.base) {
    v <- cbind(v,base)
    i.x <- rbind(c(2:n,1),1:n,rep(n+2,n))
    if (qmesh)  ## add base twice
      i.x <-  rbind(i.x,rep(n+2,n))
    i <- cbind(i,i.x)
  }
  if (qmesh) v <- rbind(v,rep(1,ncol(v))) ## homogeneous
  if (!qmesh)
    triangles3d(v[1,i],v[2,i],v[3,i],...)
  else
    return(rotate3d(qmesh3d(v,i,material=...), matrix=trans))
}

plotModel<-function(migrationIndividual, parameterVector=NULL) {
	open3d()
	if(is.null(parameterVector)) {
		parameterVector.names<-msIndividualParameters(migrationIndividual)
		parameterVector<-sequence(length(parameterVector.names)) #just dummy variables, so we can see that models are different
		names(parameterVector)<-parameterVector.names
	}
	n0multiplierParameters<-parameterVector[grep("n0multiplier",names(parameterVector))]
	migrationParameters<-parameterVector[grep("migration",names(parameterVector))]
	collapseParameters<-parameterVector[grep("collapse",names(parameterVector))]
	#now normalize for plotting
	try(n0multiplierParameters <-n0multiplierParameters/max(n0multiplierParameters))
	try(migrationParameters <-migrationParameters/max(migrationParameters))
	try(collapseParameters <-collapseParameters/max(collapseParameters))

	num.steps<-10
	base.radius<-0.3
	collapseMatrix<-migrationIndividual$collapseMatrix
	complete<-migrationIndividual$complete
	n0multiplierMap<-migrationIndividual$n0multiplierMap
	nPop<-sum(!is.na(collapseMatrix[,1]))
	current.terminal.pos<-cbind(x=1, y=0)

	for (i in sequence(nPop-1)) {
		current.terminal.pos <-rbind(current.terminal.pos, cbind(x=cos(2*pi*i/nPop), y=sin(2*pi*i/nPop)))	
	}
	text3d(x=current.terminal.pos[,1], y=current.terminal.pos[,2], z=rep(1-.2, nPop), texts=sequence(nPop))
	for (regime in sequence(dim(collapseMatrix)[2])) {
		next.terminal.pos<-current.terminal.pos*NA
		for (i in sequence(dim(collapseMatrix)[1])) {
			if (!is.na(collapseMatrix[i, regime])) {
				if(collapseMatrix[i, regime]==0) {
					next.terminal.pos[i,]<-current.terminal.pos[i,]
				}
				if(collapseMatrix[i, regime]==1) {
					next.terminal.pos[i,]<-colMeans(current.terminal.pos[which(collapseMatrix[,regime]==1),])
				}
				center<-cbind(rbind(current.terminal.pos[i,], current.terminal.pos[i,]), z=c(regime,regime+1))
				center<-cbind(seq(from=center[1,1],to=center[2,1], length.out=num.steps), seq(from=center[1,2],to=center[2,2], length.out=num.steps), seq(from=center[1,3],to=center[2,3], length.out=num.steps))
				shade3d(addNormals(subdivision3d(cylinder3d(center, radius=base.radius* n0multiplierParameters[n0multiplierMap[i, regime] ], sides=8, twist=1, closed=-2), depth=2)), col="gray")
			}
		}
		if(max(collapseMatrix[, regime],na.rm=TRUE)==1) {
			center<-cbind(rbind(current.terminal.pos[which(collapseMatrix[,regime]==1),], current.terminal.pos[which(collapseMatrix[,regime]==1),]), z=c(regime+1,regime+1))
			for (length.out in 2:6) {
			shade3d(addNormals(subdivision3d(cylinder3d(cbind(seq(from=center[1,1],to=center[2,1], length.out=length.out), seq(from=center[1,2],to=center[2,2], length.out=length.out), seq(from=center[1,3],to=center[2,3], length.out=length.out))
, radius=base.radius/3, sides=8, twist=1, closed=-2), depth=2)), col="gray") #hack to deal with rgl's twisting of cylinders
			}

		}
		rgl.viewpoint(userMatrix=matrix(c(-0.112414613366127,0.471968173980713,-0.874419331550598,0,0.993038892745972,0.0222157146781683,-0.11567322909832,0,-0.035168319940567,-0.881335437297821,-0.471180021762848,0,0,0,0,1),ncol=4))
		local.migrationMatrix<-migrationIndividual$migrationArray[,,regime]
		if(max(local.migrationMatrix, na.rm=TRUE)>0) {
			for (from.index in sequence(dim(local.migrationMatrix)[1])) {
				for(to.index in sequence(dim(local.migrationMatrix)[2])) {
					if(from.index != to.index) {
						if(!is.na(local.migrationMatrix[from.index, to.index])) {
							if(local.migrationMatrix[from.index, to.index]>0) {
								offset=1/3
								if(from.index>to.index) {
									offset=2/3	
								}
								#arrow3d(base=c(current.terminal.pos[from.index,], regime+offset), tip=c(current.terminal.pos[to.index,],regime+offset),rad=migrationParameters[local.migrationMatrix[from.index, to.index] ]*base.radius/3)
								arrow3d2(base=c(current.terminal.pos[from.index,], regime+offset), tip=c(current.terminal.pos[to.index,],regime+offset),rad=migrationParameters[local.migrationMatrix[from.index, to.index] ]*base.radius/3)
							}	
						}	
					}
				} 	
			}
		}

		current.terminal.pos<-next.terminal.pos
	}
	
}

plotModel(migrationArray[[500]], NULL)