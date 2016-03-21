###################################################################################
#
# Mark Cembrowski, April 8 2015
#
# This script looks for correlated differences in gene expression across the
# dorsal-ventral axis for the elements of the trisynaptic loop.
#
# 100415: edited to allow for count data analysis.
# 110415: edited to proofread code.
#
# numRand=10: numeric.  Integer number of stochastic trials.
#
###################################################################################

correlatedDv <- function(fpkmThres=10,foldThres=2,numRand=10,quiet=T,export=F,
				countThres=20,doFpkm=T){
	if(doFpkm){
		theThres <- fpkmThres
	}else{
		theThres <- countThres
	}


	# Extract the lists that correspond to the dorsal and ventral enriched genes.
	toMask <- colnames(fpkmPoolMat)[!grepl('ca1_d|ca1_v',colnames(fpkmPoolMat))]
	ca1dE <- enrichedGenes('ca1_d',fpkmThres=theThres,foldThres=foldThres,
			avgPass=T,mask=toMask,doFpkm=doFpkm)
	ca1vE <- enrichedGenes('ca1_v',fpkmThres=theThres,foldThres=foldThres,
			avgPass=T,mask=toMask,doFpkm=doFpkm)

	toMask <- colnames(fpkmPoolMat)[!grepl('ca3_d|ca3_v',colnames(fpkmPoolMat))]
	ca3dE <- enrichedGenes('ca3_d',fpkmThres=theThres,foldThres=foldThres,
			avgPass=T,mask=toMask,doFpkm=doFpkm)
	ca3vE <- enrichedGenes('ca3_v',fpkmThres=theThres,foldThres=foldThres,
			avgPass=T,mask=toMask,doFpkm=doFpkm)

	toMask <- colnames(fpkmPoolMat)[!grepl('dg_d|dg_v',colnames(fpkmPoolMat))]
	dgdE <- enrichedGenes('dg_d',fpkmThres=theThres,foldThres=foldThres,
			avgPass=T,mask=toMask,doFpkm=doFpkm)
	dgvE <- enrichedGenes('dg_v',fpkmThres=theThres,foldThres=foldThres,
			avgPass=T,mask=toMask,doFpkm=doFpkm)
	
	if(doFpkm){
		theMat <- fpkmPoolMat
	}else{
		theMat <- cpmPoolMat
	}
	
	# Get total genes expressed in each group.
	ca1d <- rownames(theMat)[theMat$ca1_d>theThres]
	ca1v <- rownames(theMat)[theMat$ca1_v>theThres]

	ca3d <- rownames(theMat)[theMat$ca3_d>theThres]
	ca3v <- rownames(theMat)[theMat$ca3_v>theThres]

	dgd <- rownames(theMat)[theMat$dg_d>theThres]
	dgv <- rownames(theMat)[theMat$dg_v>theThres]

	# Get number enriched in >=2 comparisons.
	ca1ca3dE <- intersect(ca1dE,ca3dE)
	ca1ca3vE <- intersect(ca1vE,ca3vE)
	ca1dgdE <- intersect(ca1dE,dgdE)
	ca1dgvE <- intersect(ca1vE,dgvE)
	ca3dgdE <- intersect(ca3dE,dgdE)
	ca3dgvE <- intersect(ca3vE,dgvE)
		
	# Get number enriched in 3 comparisons.
	ca1ca3dgdE <- intersect(intersect(ca1dE,ca3dE),dgdE)
	ca1ca3dgvE <- intersect(intersect(ca1vE,ca3vE),dgvE)

	# Convert gene lists to numbers.
	nca1d <- length(ca1d)
	nca3d <- length(ca3d)
	ndgd <- length(dgd)
	nca1v <- length(ca1v)
	nca3v <- length(ca3v)
	ndgv <- length(dgv)

	print('The number of genes expressed at each location: DGd/CA3d/CA1d DGv/CA3v/CA1v')
	print(paste(ndgd,nca3d,nca1d,ndgv,nca3v,nca1v))	

	nca1dE <- length(ca1dE)
	nca3dE <- length(ca3dE)
	ndgdE <- length(dgdE)
	nca1vE <- length(ca1vE)
	nca3vE <- length(ca3vE)
	ndgvE <- length(dgvE)
	
	print('The number of genes enrihed in each pairwise comparison: DGd/CA3d/CA1d DGv/CA3v/CA1v')
	print(paste(ndgdE,nca3dE,nca1dE,ndgvE,nca3vE,nca1vE))
	
	nca1ca3dE <- length(ca1ca3dE)
	nca1dgdE <- length(ca1dgdE)
	nca3dgdE <- length(ca3dgdE)
	nca1ca3vE <- length(ca1ca3vE)
	nca1dgvE <- length(ca1dgvE)
	nca3dgvE <- length(ca3dgvE)

	nca1ca3dgdE <- length(ca1ca3dgdE)
	nca1ca3dgvE <- length(ca1ca3dgvE)
	
	# Draw some guys by chance.
	theOut <- as.data.frame(matrix(nrow=numRand,ncol=8))
	colnames(theOut) <- c('ca1ca3d','ca1dgd','ca3dgd','ca1ca3dgd',
			      'ca1ca3v','ca1dgv','ca3dgv','ca1ca3dgv')

	for (ii in 1:numRand){
		theOut[ii,1] <- .randomOverlap(ca1d,nca1dE,ca3d,nca3dE)
		theOut[ii,2] <- .randomOverlap(ca1d,nca1dE,dgd,ndgdE)
		theOut[ii,3] <- .randomOverlap(ca3d,nca3dE,dgd,ndgdE)
		theOut[ii,4] <- .randomOverlap(ca1d,nca1dE,ca3d,nca3dE,dgd,ndgdE)
		
		theOut[ii,5] <- .randomOverlap(ca1v,nca1vE,ca3v,nca3vE)
		theOut[ii,6] <- .randomOverlap(ca1v,nca1vE,dgv,ndgvE)
		theOut[ii,7] <- .randomOverlap(ca3v,nca3vE,dgv,ndgvE)
		theOut[ii,8] <- .randomOverlap(ca1v,nca1vE,ca3v,nca3vE,dgv,ndgvE)
	}
		
	muOut <- apply(theOut,2,mean)
	sdOut <- apply(theOut,2,sd)
	errOut <- sdOut*2
	errLo <- muOut-errOut
	errHi <- muOut+errOut

	expData <- as.data.frame(
		c(nca1ca3dE,nca1dgdE,nca3dgdE,nca1ca3dgdE,
		  nca1ca3vE,nca1dgvE,nca3dgvE,nca1ca3dgvE))
	colnames(expData) <- 'expData'
	rownames(expData) <- colnames(theOut)
	
	if(!quiet){
		print(expData)
		for (ii in 1:8){
			# Print out p values for everything.
			print(sum(theOut[,ii]>expData[ii,1])/numRand)
		}
	}

return(theOut)

	df <- as.data.frame(cbind(rownames(expData),muOut,errLo,errHi,expData))
	colnames(df)[1] <- 'test'

	df$test <- factor(df$test,levels=c('ca3dgd','ca1ca3d','ca1dgd','ca1ca3dgd',
		'ca3dgv','ca1ca3v','ca1dgv','ca1ca3dgv'))
	
	gg <- ggplot(df,aes(x=test)) + geom_point(aes(y=expData),colour='red',size=12,shape=95)
	gg <- gg + geom_errorbar(aes(ymin=errLo,ymax=errHi),width=0.5)
	gg <- gg + geom_point(aes(y=muOut),colour='black',size=4,shape=95)
	gg <- gg + theme_bw()
	print(gg)

	if(!quiet){
		print('here are mean and standard deviation')
		print(muOut)
		print(sdOut)
	}

	if(!quiet){
		print('These are the genes enriched dorsally across all comparisons:')
		print(as.data.frame(sort(idToSym(ca1ca3dgdE))))
		print('These are the genes enriched ventrally across all comparisons:')
		print(as.data.frame(sort(idToSym(ca1ca3dgvE))))
	}

	theD <- theMat[ca1ca3dgdE,]
	theD <- cbind(theD,rep('dorsal',nrow(theD)))
	colnames(theD)[ncol(theD)] <- 'regionEnriched'
	theD <- theD[,c('regionEnriched','dg_d','ca3_d','ca1_d','dg_v','ca3_v','ca1_v')]
	rownames(theD) <- idToSym(rownames(theD))
	theD <- theD[order(rownames(theD)),]
	
	theV <- theMat[ca1ca3dgvE,]
	theV <- cbind(theV,rep('ventral',nrow(theV)))
	colnames(theV)[ncol(theV)] <- 'regionEnriched'
	theV <- theV[,colnames(theD)]
	rownames(theV) <- idToSym(rownames(theV))
	theV <- theV[order(rownames(theV)),]

	theOut <- rbind(theD,theV)

	# Export data for checking into DAVID, if desired.
	if(export){
		print('here')
		write.table(c(ca1ca3dgdE,ca1ca3dgvE),file='dvCorrGenes.txt',sep='\t',quote=F,row.names=F)
	}


	invisible(theOut)

#	print('Number of genes expressed in CA1 / enriched in CA1 (dorsal vs ventral):')
#	print(paste(length(ca1d),length(ca1dE),length(ca1v),length(ca1vE)))
#	print('Number of genes expressed in CA3 / enriched in CA3 (dorsal vs ventral):')
#	print(paste(length(ca3d),length(ca3dE),length(ca3v),length(ca3vE)))
#	print('Number of genes expressed in DG / enriched in DG (dorsal vs ventral):')
#	print(paste(length(dgd),length(dgdE),length(dgv),length(dgvE)))

}

# Draw numA genes from list geneA and numB genes from list geneB and 
# calculate the number of intersections.  Optional arguments provide
# a third list to supply.
.randomOverlap <- function(geneA,numA,geneB,numB,geneC=c(),numC=-1){

	fromA <- sample(geneA,numA)
	fromB <- sample(geneB,numB)
	if(numC>0.1){
		fromC <- sample(geneC,numC)
		invisible(length(intersect(intersect(fromA,fromB),fromC)))
	}else{
		invisible(length(intersect(fromA,fromB)))
	}
}
