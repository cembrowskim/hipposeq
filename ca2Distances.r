#######################################################################################
#
# Mark Cembrowski, Janelia Research Campus, Jan 25 2016
#
# This script looks for the Euclidean distance between the CA1, CA2, and CA3 datasets,
# as requested by Reviewer 2.
#
#######################################################################################

ca2Distances <- function(doLog=T,pseudo=1){

	cadMat <- grepl('ca3_d|ca2|ca1_d',colnames(fpkmRepMat))
	cadMat <- fpkmRepMat[,cadMat]
	cadMat <- t(cadMat)

	if(doLog){
		cadMat <- cadMat+pseudo
		cadMat <- log10(cadMat)
	}

	theDist <- as.matrix(dist(cadMat,method='euclidean',diag=T,upper=T))

	ca2ca3 <- theDist[4:6,1:3]
	ca2ca1 <- theDist[4:6,7:9]

	ca2ca1mean <- mean(ca2ca1)
	ca2ca1err <- sd(ca2ca1)

	ca2ca3mean <- mean(ca2ca3)
	ca2ca3err <- sd(ca2ca3)

	print(paste('CA2-CA1:',ca2ca1mean,'+/-',ca2ca1err))
	print(paste('CA2-CA3:',ca2ca3mean,'+/-',ca2ca3err)) 

	print(t.test(ca2ca3,ca2ca1,paired=T))
}
