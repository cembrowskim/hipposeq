######################################################################################################
# Mark Cembrowski, Janelia Farm, Sept 13 2015
#
# This script plots the number of enriched genes as a function of fold change for two datasets.
#
# 100415: edited to work for count data.
#
# OPTIONAL INPUTS: 
# fpkmThres=10: variable.  Threshold for cutting off genes.
# rMax=-1: variable.  Maximum y value on rendered graph.
#
######################################################################################################

enrichmentVsFold <- function(popA,popB,fpkmThres=10,cpmThres=20,doFpkm=T,rMax=-1){
	# Retrieve all values.
	if(doFpkm){
		theMat <- fpkmPoolMat
		theThres <- fpkmThres
	}else{
		theMat <- cpmPoolMat
		theThres <- cpmThres
	}

	aEnr <- .getFoldDiff(theMat,popA,popB,theThres)
	bEnr <- .getFoldDiff(theMat,popB,popA,theThres)
	abEnr <- aEnr
	abEnr$yAxis <- abEnr$yAxis + bEnr$yAxis
	abEnr$yAxisNorm <- abEnr$yAxis/(abEnr$yAxis[1])*100
	
	# Do a plot for two lines for each comparison (e.g., one line for dosally-enriched transcripts
	# and one line for ventrally-enriched transcripts).
	plotOut <- ggplot(aEnr,aes_string(x='xAxis',y='yAxis'))
	plotOut <- plotOut + geom_line(colour='green')
	plotOut <- plotOut + geom_line(data=bEnr,aes_string(x='xAxis',y='yAxis'),colour='purple')
	plotOut <- plotOut + geom_line(data=abEnr,aes_string(x='xAxis',y='yAxis'),colour='black')
	plotOut <- plotOut + labs(x='Fold change',y='Number of genes')
	if(rMax<0){
		plotOut <- plotOut + ylim(c(0,abEnr$yAxis[1]))
	}else{
		plotOut <- plotOut + ylim(c(0,rMax))
	}
	plotOut <- plotOut + theme_bw()
	print(plotOut)


}

# Single-serving function filtering the original FPKM matrix.
.getFoldDiff <- function(fpkmMat,sampleX,sampleY,fpkmThres){
	procMat <- fpkmMat[,c(sampleX,sampleY)]
	toFilt <- procMat[,1]>fpkmThres
	print(paste(sum(toFilt),'of',nrow(procMat),'passed threshold'))

	procMat <- procMat[toFilt,]

	foldMat <- procMat[,1] / procMat[,2]

	xHigh <- 10
	xAxis <- 3:xHigh
	yAxis <- xAxis*0
	for (ii in 1:length(xAxis)){
		curX <- xAxis[ii]
		curY <- round(sum(foldMat>curX))
		yAxis[ii] <- curY
	}
	yAxisNorm <- yAxis/yAxis[1]*100

	outList <- data.frame(xAxis,yAxis,yAxisNorm)

	return(outList)
}
