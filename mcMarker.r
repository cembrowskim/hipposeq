###########################################################################################
#
# Mark Cembrowski, Janelia Research Campus, Apr 7 2015
# 
# This script plots mossy cell marker genes.
#
###########################################################################################

mcMarker <- function(doFpkm=T,oldMarks=T,quiet=F){
	if(oldMarks){
		mcMarksOld <- c('Calb2') # need to look up old references on this ...
	}else{
		mcMarksOld <- c()
	}

	if(!doFpkm){
		theThres <- 0.01 # Small, just validating >3 fold different
	}else{
		theThres <- 10
	}

	mcMarksNew <- enrichedGenes('ca4',avgPass=T,fpkmThres=theThres,foldThres=3,doFpkm=doFpkm)
	mcMarksNew <- idToSym(mcMarksNew)
	mcMarksNew <- mcMarksNew[order(mcMarksNew)]

	if(!quiet){	
		print('These are the genes to be looked at: ')
		print(mcMarksNew)
	}
	
	mcMarks <- c(mcMarksOld,mcMarksNew)
	if(!quiet){
		print(paste('The number of new marker genes is ',length(setdiff(mcMarks,mcMarksOld))))
		print(fpkmHeatmap(mcMarks,doNorm=T,replicates=T))
	}

	invisible(mcMarks)
}

