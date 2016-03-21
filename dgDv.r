#############################################################################
#
# Mark Cembrowski, Janelia Research Campus, April 3 2015
#
# This script plots the number of D-V enriched genes in DG granule cells.
#
#############################################################################

dgDv <- function(fpkmThres=10){

	theFolds <- c(3:10)
	df <- as.data.frame(matrix(nrow=theFolds,ncol=2))
	for (ii in 1:length(theFolds)){
		df[ii,] <- .dgFold(theFolds[ii],fpkmThres)
	}
	total <- df[,1]+df[,2]
	df <- cbind(theFolds,df,total)
	rownames(df) <- theFolds
	colnames(df) <- c('theFolds','dorsal','ventral','total')

	print(df)

	# Plot.
	gg <- ggplot(df,aes(x=theFolds))
	gg <- gg + geom_line(aes(y=dorsal),colour='green')
	gg <- gg + geom_line(aes(y=ventral),colour='magenta')
	gg <- gg + geom_line(aes(y=total),colour='black')
	gg <- gg + theme_bw()
	gg <- gg + expand_limits(y=0)
	gg <- gg + xlab('Fold difference') + ylab('Number of genes')
	print(gg)	
}

# details = F: return only numbers of genes, or gene names as well?
.dgFold <- function(theFold,fpkmThres,details=F){
	# Get genes that are enriched in either direction with the fold thres.
	dgMat <- fpkmPoolMat[,c('dg_d','dg_v')]
	dgThres <- apply(dgMat,1,max)>fpkmThres
	dgMat <- dgMat[dgThres,]	

	dgRat <- dgMat$dg_d/dgMat$dg_v
	

	dgD <- sum(dgRat>theFold)
	dgDNames <- idToSym(rownames(dgMat)[dgRat>theFold])
	dgV <- sum(dgRat<(1/theFold))
	dgVNames <- idToSym(rownames(dgMat)[dgRat<(1/theFold)])


	if(!details){
		invisible(c(dgD,dgV))
	}else{
		dgNames <- c(dgDNames,dgVNames)
		dgEnriched <- c(rep('dorsal',length(dgDNames)),rep('ventral',length(dgVNames)))
		theOut <- cbind(dgNames,dgEnriched)

		dgFpkm <- fpkmPoolMat[rownames(theOut),c('dg_d','dg_v')]
		dgFold <- dgFpkm$dg_d/dgFpkm$dg_v
		dgFold <- pmax(dgFold,1/dgFold)

		theOut <- cbind(theOut,dgFpkm,dgFold)

		# Put into order of fold change.
		theOutD <- subset(theOut,dgEnriched=='dorsal')
		theOutD <- theOutD[order(-theOutD$dgFold),]
		theOutV <- subset(theOut,dgEnriched=='ventral')
		theOutV <- theOutV[order(-theOutV$dgFold),]
		theOut <- rbind(theOutD,theOutV)

		invisible(theOut)
	}
}	
