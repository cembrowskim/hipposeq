#################################################################################
#
# Mark Cembrowski, Janelia Research Campus, Apr 7 2015
#
# This looks for the number of genes enriched in CA2 relative to the rest of the
# principal cells of the hippocampus.
#
#################################################################################

ca2EnrichVsFold <- function(fpkmThres=10){

	theFolds <- c(2:10)
	df <- as.data.frame(matrix(nrow=theFolds,ncol=1))
	for (ii in 1:length(theFolds)){
		numCa2Enrich <- length(enrichedGenes('ca2',fpkmThres=fpkmThres,
			foldThres=theFolds[ii],avgPass=T))
		df[ii,] <- numCa2Enrich
	}
	
	df <- cbind(theFolds,df)
	colnames(df) <- c('theFolds','total')

	print(df)

	# Plot.
	gg <- ggplot(df,aes(x=theFolds))
	gg <- gg + geom_line(aes(y=total),colour='black')
	gg <- gg + theme_bw()
	gg <- gg + expand_limits(y=0)
	gg <- gg + xlab('Fold difference') + ylab('Number of genes')
	print(gg)	

	invisible(gg)
}
