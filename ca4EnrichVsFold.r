####################################################################################
#
# Mark Cembrowski, Janelia Research Campus, April 8 2015
#
# This script looks for the number of mossy cell-enriched transcripts relative
# to neighbouring DG GCs and CA3 PCs.
#
####################################################################################

ca4EnrichVsFold <- function(fpkmThres=10){

	theFolds <- c(2:10)
	df <- as.data.frame(matrix(nrow=theFolds,ncol=3))

	theMask <- setdiff(colnames(fpkmPoolMat),c('dg_d','ca4','ca3_d'))

	for (ii in 1:length(theFolds)){
		numDgEnrich <- length(enrichedGenes('dg_d',fpkmThres=fpkmThres,
			foldThres=theFolds[ii],avgPass=T,mask=theMask))
		numCa4Enrich <- length(enrichedGenes('ca4',fpkmThres=fpkmThres,
			foldThres=theFolds[ii],avgPass=T,mask=theMask))
		numCa3Enrich <- length(enrichedGenes('ca3_d',fpkmThres=fpkmThres,
			foldThres=theFolds[ii],avgPass=T,mask=theMask))
		df[ii,] <- c(numDgEnrich,numCa4Enrich,numCa3Enrich)
	}
	
	df <- cbind(theFolds,df)
	colnames(df) <- c('theFolds','dg','ca4','ca3')

	# Plot.
	gg <- ggplot(df,aes(x=theFolds))
	gg <- gg + geom_line(aes(y=dg),colour='red')
	gg <- gg + geom_line(aes(y=ca4),colour='magenta')
	gg <- gg + geom_line(aes(y=ca3),colour='green')
	gg <- gg + theme_bw()
	gg <- gg + expand_limits(y=0)
	gg <- gg + xlab('Fold difference') + ylab('Number of genes')
	print(gg)	

	invisible(df)
}
