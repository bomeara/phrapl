library(rgl)
library(phrapl)
library(RColorBrewer)
load("/Users/bomeara/Documents/MyDocuments/Active/phrapl/pkg/data/migrationArray_npop3_maxK5.Rsave")


arrow3d2 <- function(base, tip, rad=1, col="red",offset=TRUE) {
	length.out=10
	if(!offset) {
		length.out=8	
	}
	new.pos<-cbind(seq(from=base[1], to=tip[1], length.out=length.out), seq(from=base[2], to=tip[2], length.out=length.out), seq(from=base[3], to=tip[3], length.out=length.out))
                   	#base<-new.pos[6,]
                   	#tip<-new.pos[8,]
                   #	lines3d(x=new.pos[1:8,1], y=new.pos[1:8,2], z=new.pos[1:8,3], col=col, lwd=rad*20)

		shade3d(addNormals(subdivision3d(cylinder3d(rbind(new.pos[8,], new.pos[6,]), radius=c(0.00001, rad), sides=8, twist=1, closed=-2), depth=2)), col=col)
	shade3d(addNormals(subdivision3d(cylinder3d(rbind(new.pos[7,], new.pos[1,]), radius=rep(rad/4,2), sides=8, twist=1, closed=-2), depth=2)), col=col)
}

#suggestion from Zach Marion: allow placement on a landscape map, so you can see rivers and such as boundaries
plotModel<-function(migrationIndividual, parameterVector=NULL, taxonNames=NULL, time.axis=FALSE, time.axis.col="black", apply.base=FALSE, base.color='black', new.window=TRUE) {
if(new.window) {
	open3d()
	}
  rgl.material(shininess=100)
	
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
	if(is.null(taxonNames)) {
	  taxonNames<-sequence(nPop)
	}
	
	for (i in sequence(nPop-1)) {
		current.terminal.pos <-rbind(current.terminal.pos, cbind(x=cos(2*pi*i/nPop), y=sin(2*pi*i/nPop)))	
	}
	text3d(x=current.terminal.pos[,1], y=current.terminal.pos[,2], z=rep(1-.2, nPop), texts=taxonNames, col="black")
	if(time.axis) {
		arrow3d2(base=c(2,2,1+dim(collapseMatrix)[2]), tip=c(2,2,1), rad=0.1, col=time.axis.col, offset=FALSE)
		text3d(x=2, y=2, z=.8, texts="Present day", col=time.axis.col, cex=0.8)
	}
	if (apply.base) {
		plane.coords<-combn(apply(expand.grid(x=c(-2,2), y=c(-2,2), z=c((1+dim(collapseMatrix)[2]),(1.25+dim(collapseMatrix)[2]))), 1, paste, collapse="_"),4)
		for (plane.element in sequence(dim(plane.coords)[2])) {
			corners<-matrix(as.numeric(unlist(strsplit(plane.coords[,plane.element], "_"))), byrow=TRUE, ncol=3)
			quads3d(x=corners, col=base.color)
		}
	}
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
				center<-cbind(rbind(current.terminal.pos[i,], current.terminal.pos[i,]), z=c(regime,regime+1.2))
				center<-cbind(seq(from=center[1,1],to=center[2,1], length.out=num.steps), seq(from=center[1,2],to=center[2,2], length.out=num.steps), seq(from=center[1,3],to=center[2,3], length.out=num.steps))
				shade3d(addNormals(subdivision3d(cylinder3d(center, radius=base.radius* n0multiplierParameters[n0multiplierMap[i, regime] ], sides=8, twist=1, closed=-2), depth=2)), col=c("gray",rev(brewer.pal(8,"Set1")))[n0multiplierMap[i, regime]])
			}
		}
		if(max(collapseMatrix[, regime],na.rm=TRUE)==1) {
			center<-cbind(rbind(current.terminal.pos[which(collapseMatrix[,regime]==1),], current.terminal.pos[which(collapseMatrix[,regime]==1),]), z=c(regime+1,regime+1.2))
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
								arrow3d2(base=c(current.terminal.pos[from.index,], regime+offset), tip=c(current.terminal.pos[to.index,],regime+offset),rad=migrationParameters[local.migrationMatrix[from.index, to.index] ]*base.radius/3, col=(brewer.pal(8,"Set1"))[local.migrationMatrix[from.index, to.index]])
							}	
						}	
					}
				} 	
			}
		}

		current.terminal.pos<-next.terminal.pos
	}
	
}

saveMovie<-function(total.revolutions=1, duration=10, save.dir=NULL) {
	dir=tempdir()
	if(!is.null(save.dir)) {
		try(system(paste("mkdir", save.dir)), silent=TRUE)	
		dir=save.dir
	}
	result<-movie3d(spin3d(rpm=total.revolutions/(duration/60)), duration=duration, dir=dir)
}

plotModel(migrationArray[[500]], taxonNames=c("","",""))
#saveMovie(save.dir="~/Desktop/movies")

try(system("mkdir ~/Desktop/models"))
setwd("~/Desktop/models")
for (i in sequence(length(migrationArray))) {
	plotModel(migrationArray[[i]], taxonNames=c("","",""), new.window=FALSE)
	rgl.snapshot(paste("model",i,".png", sep=""))	
	saveMovie(save.dir="~/Desktop/models")
	system(paste("mv movie.gif movie", i, ".gif", sep=""))
	rgl.clear()
}