#################################################################################
#
# Mark Cembrowski, Janelia Research Campus, Apr 7 2015
#
# This looks for the number of genes enriched in the various CA region sof the 
# hippocampus.
#
#################################################################################

caRegionEnrichVsFold <- function(fpkmThres=10,countThres=20,doFpkm=T){
	if(doFpkm){
		theThres <- fpkmThres
	}else{
		theThres <- countThres
	}

	theFolds <- c(3:10)
	df <- as.data.frame(matrix(nrow=theFolds,ncol=3))

	theMask <- setdiff(colnames(fpkmPoolMat),c('ca3_d','ca2','ca1_d'))

	for (ii in 1:length(theFolds)){
		numCa3Enrich <- length(enrichedGenes('ca3_d',fpkmThres=theThres,
			foldThres=theFolds[ii],avgPass=T,mask=theMask,doFpkm=doFpkm))
		numCa2Enrich <- length(enrichedGenes('ca2',fpkmThres=theThres,
			foldThres=theFolds[ii],avgPass=T,mask=theMask,doFpkm=doFpkm))
		numCa1Enrich <- length(enrichedGenes('ca1_d',fpkmThres=theThres,
			foldThres=theFolds[ii],avgPass=T,mask=theMask,doFpkm=doFpkm))
	if(theFolds[ii]>9.5){
print(idToSym(enrichedGenes('ca1_d',fpkmThres=theThres,
                        foldThres=theFolds[ii],avgPass=T,mask=theMask,doFpkm=doFpkm)))
}
		df[ii,] <- c(numCa3Enrich,numCa2Enrich,numCa1Enrich)
	}
	
	df <- cbind(theFolds,df)
	colnames(df) <- c('theFolds','ca3','ca2','ca1')

	# Plot.
	gg <- ggplot(df,aes(x=theFolds))
	gg <- gg + geom_line(aes(y=ca2),colour='black')
	gg <- gg + geom_line(aes(y=ca3),colour='red')
	gg <- gg + geom_line(aes(y=ca1),colour='blue')
	gg <- gg + theme_bw()
	gg <- gg + expand_limits(y=0)
	gg <- gg + xlab('Fold difference') + ylab('Number of genes')
	gg <- gg + ggtitle('Blue:CA1,Black:CA2,red:CA3')
	print(gg)	

	invisible(df)
}
