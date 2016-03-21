###################################################################################
#
# Mark Cembrowski, Janelia Farm, October 10 2013
#
# This script implements a whole genome coexpression network analysis (WGCNA)
# as discussed in Zhang and Horvath 2005.  The top script is wgcna, and calls the
# remainder of the scripts.
#
# INPUT (for wgcna):
# numGenes: integer.  Number of genes to run the analysis on [selected for
#	according to their variance across the dataset] 
# cutHeight: numeric.  Height to cut tree to form clusters.  If <0, plots the
#	dendrogram and then halts execution, allowing the user to visualise
#	where they want the cut to be.
# 
# OUTPUT:
# 1. The cut groups as saved, in the directory 'wgcna/cutXgenesX/'
# 2. A plot is rendered (see optional inputs).
#
###################################################################################

# If cut height <0, execution is halted after dendrogram generated, allowing user
# to set cut height.
wgcna <- function(numGenes,cutHeight,thePlot='dend'){
	corOut <- .getCorMat(numGenes) # Get correlation matrix to use; 'S' in Z&H.
	disOut <- .corToDis(corOut) # Get dissimilarity measure.
	clusOut <- .groupDis(disOut,cutHeight) # Generate clusters.
	.plotFeatures(disOut,clusOut) # Visualise.
	.saveClusters(clusOut,cutHeight,numGenes)
}

.saveClusters <- function(clusts,cutHeight,numGenes){
	theClusts <- unique(clusts$groups)
	theClusts <- theClusts[theClusts>0] # Remove nonclustered elements.

	# Remove any existing text files in the next save directory; this is to
	# prevent previous files from accidentally sticking around on a new run.
	saveDir <- paste('cut',cutHeight,'genes',numGenes,sep='')
	dir.create(file.path(getwd(),'wgcna'),showWarnings=F) # Create wgcna dir
	dir.create(file.path(getwd(),'wgcna',saveDir),showWarnings=F) # Create directory.

	txtFiles <- list.files(saveDir,pattern='.txt')
	for (txtFile in txtFiles){
		file.remove(paste(saveDir,txtFile,sep=''))
	}

	for (ii in theClusts){
		geneCohort <- names(clusts$groups[clusts$groups==ii])
		getEnsId(geneCohort,saveName=paste('clust',ii,sep=''),
			saveDir=paste('wgcna',saveDir,sep='/'))
	}

	# David multilist that shit.
	.davidMultilist(paste('wgcna',saveDir,sep='/'))
}


# Plot features of the WGNCA; namely, the TOM, groupings, and heat maps.
.plotFeatures <- function(theDis,clusts){
	# Clean up dissimilarity matrix for plotting.
	sorted <- theDis[,clusts$dend$order.lab]
	sorted <- sorted[clusts$dend$order.lab,]
	meltSorted <- melt(sorted)
	meltSorted$Var1 <- factor(meltSorted$Var1,levels=clusts$dend$order.lab)
	meltSorted$Var2 <- factor(meltSorted$Var2,levels=clusts$dend$order.lab)

	# Plot dissimilarity matrix.
	corPlot <- ggplot(meltSorted,aes_string(x='Var1',y='Var2'))
	corPlot <- corPlot + geom_tile(aes_string(fill='value'))
	corPlot <- corPlot + scale_fill_continuous(low="red",high="white")
	corPlot <- corPlot + theme(text = element_blank())
	corPlot <- corPlot + theme(axis.ticks = element_blank())
	corPlot <- corPlot + theme(legend.position='none')

	# Prep a data frame for showing data frame enrollment.
	groupLine <- data.frame(gene_id=names(clusts$groups),grouping=clusts$groups,preFilt=clusts$groupsPre)
	groupLine$gene_id <- factor(groupLine$gene_id,levels=clusts$dend$order.lab)

	# Reorganise the grouping numbers to enhance contrast between neighbouring groups
	# in plots; note that this does not change the group values anywhere else, but
	# is used strictly for visualisation.
	temper <- as.vector(groupLine$preFilt)
	temper[temper%%2<0.1] <- temper[temper%%2<0.1] + max(temper)
	groupLine$preFilt <- temper
	temper <- as.vector(groupLine$grouping)
	temper[temper%%2<0.1] <- temper[temper%%2<0.1] + max(temper)
	groupLine$grouping <- temper

	# Plot a colourised line that corresponds to group enrollment.
	groupPlot <- ggplot(groupLine,aes_string(x='gene_id',y='1')) + geom_tile(aes_string(fill='factor(grouping)'))
	groupPlot <- groupPlot + theme(axis.title.x=element_blank(),axis.text = element_blank())
	groupPlot <- groupPlot + labs(x='Gene',y='Group')
	groupPlot <- groupPlot + theme(axis.ticks = element_blank())
	groupPlot <- groupPlot + theme(legend.position='none')
	groupPlot <- groupPlot + theme(panel.grid=element_blank())
	groupPlot <- groupPlot + theme(panel.background=element_blank())

	# Plot a similar line that corresponds to group enrollment, before applying
	# minimum number of entries.
	groupPrePlot <- ggplot(groupLine,aes_string(x='gene_id',y='1')) + geom_tile(aes_string(fill='factor(preFilt)'))
        groupPrePlot <- groupPrePlot + theme(axis.title.x=element_blank(),axis.text = element_blank())
        groupPrePlot <- groupPrePlot + labs(x='Gene',y='Group') 
        groupPrePlot <- groupPrePlot + theme(axis.ticks = element_blank())
        groupPrePlot <- groupPrePlot + theme(legend.position='none')
        groupPrePlot <- groupPrePlot + theme(panel.grid=element_blank())
        groupPrePlot <- groupPrePlot + theme(panel.background=element_blank())

	# Display a heat plot of genes and the associated group enrollment.	
	heatPlot <- fpkmHeatmap(groupLine$gene_id,doNorm=T,legend.position='none',
			axis.text.y = element_blank(),axis.ticks.y = element_blank())
	par(ask=T)
	print(plot(as.dendrogram(clusts$dend)))
	print(multiplot(groupPrePlot,groupPlot,heatPlot,cols=1))

	# Do plots.
	print(corPlot)
	par(ask=F)
	
	invisible()
}

# Transform distance matrix to groups of genes by cutting the tree.
# OPTIONAL ARGUMENTS:
# numEnts=10: minimum number of entries to consider a group.
.groupDis <- function(theDis,cutHeight,numEnts=10){
	library(cluster)
	theDend <- diana(theDis)
	
	# If cutHeight<0, execution is aborted; this is used for visualising
	# the dendrogram to choose a cut height.
	if(cutHeight<0){
		print(plot(as.dendrogram(theDend)))
		stop('Cut height <0, stopping b/c dendrogram plotted')
	}

	# Cut tree at a given height.
	preFilt <- cutree(as.hclust(theDend),h=cutHeight)
	preFilt <- preFilt[theDend$order.lab]

	# Relabel groups in a monotonically increasing order, makes groups
	# easier to identify in graphs.
	tempFilt <- preFilt
	origOrd <- unique(preFilt)
	newOrd <- 1:length(origOrd)
	for (ii in newOrd){
		preFilt[tempFilt==origOrd[ii]] <- ii
	}

	theCut <- preFilt

	# Obtain groups with more than numEnts elements.
	groupsPastCutoff <- which(table(theCut)>numEnts)
	
	# Map groups not making cutoff to -1.
	theCut[!theCut%in%groupsPastCutoff] <- -1

	# Order theCut in the order it would be plotted.
	theCut <- theCut[theDend$order.lab]
	preFilt <- preFilt[theDend$order.lab]

	# Return dendrogram and groupings.	
	theOut <- list(dend=theDend,groups=theCut,groupsPre=preFilt)
	return(theOut)
}
	
# Convert the correlation matrix to the dissimilarity matrix, according
# to the methods of Zhang and Horvath. For now, 1-(corr coeff)
.corToDis <- function(corMat){
	D <- 1 - corMat
	return(D)
}

# Builds original correlation matrix and transforms to a nonnegative matrix.
.getCorMat <- function(numGenes,fpkmThres=5){
	# Use FPKM cutoff matrix to screen out low-expressing genes.  Typically
	# this does not make a big difference, as the large variance genes
	# tend to have large FPKM values anyhow...
	filtMat <- fpkmMatrix(genes(cuff_data))
	filtMat <- apply(filtMat,1,max)>fpkmThres

	# Get FPKM matrix.
	theMat <- repFpkmMatrix(genes(cuff_data))
	theMat <- theMat[filtMat,]

	# Changes the FPKM matrix to a correlation coefficient matrix
	# for the most variable genes.
        theVar <- apply(theMat,1,var)
        theOrd <- order(theVar)
        toKeep <- tail(theOrd,numGenes)
        temp <- theMat[toKeep,]
        temp <- as.matrix(temp)
	
	# Make coeff between zero and one via a linear shift
	theCorr <- (1+cor(t(temp)))/2

        return(theCorr)
}

# This function plots the MDS, but is now deprecated.  MDS plots look crummy.
.plotMds <- function(groupings){
	subMat <- t(subFpkmMatrix(names(groupings)))

	d <- JSdist(makeprobs(subMat))
	fit <- cmdscale(d,eig=T,k=2)

	res <- data.frame(names=rownames(fit$points),M1=fit$points[,1],
		M2=fit$points[,2],groupings=groupings)
	p <- ggplot(res)
	p <- p + geom_point(aes(x=M1,y=M2,color=factor(groupings),size=(groupings>0+1))) + theme_bw()
	print(p)
}

