################################################################################################
#
# Mark Cembrowski, Janelia Research Campus, Oct 15 2015
#
# This function takes a list of genes and renders a heatmap.  101515: extended to enable
# CPM-based data to be visualised.
#
# INPUT:
# geneList: character.  One or more gene names to plot; can be ENS IDs or gene short names.
#
# OPTIONAL ARGUMENTS:
# doNorm=F: logical. normalise each gene to its peak value
# replicates=F: logical. show individual replicates
# mask=c(): character.  mask any supplied samples from analysis and visualisation.
#
# OUTPUT:
# A heatmap is rendered for the input genes and returned.
#
################################################################################################

geneHeatmap <- function(geneList,doNorm=F,replicates=F,mask=c(),rangeMax=1e9,doFpkm=T,...){
	# Check for duplicate entries in gene list; eliminate them if present and
	# alert user.
	if(sum(duplicated(geneList))>0.1){
		warning('Duplicated entries in gene list; removing.')
		geneList <- geneList[!duplicated(geneList)]
	}

	# Check to make sure that geneList is character vector; unintended behaviour can 
	# emerge if using factors.  If this happens, use the as.vector() wrapper.
	if(class(geneList)!='character'){
		geneList <- as.vector(geneList)
	}

	# Convert gene list to unique Ensembl IDs, if supplied with gene short names.
	if(!isGeneId(geneList)){
		geneList <- symToId(geneList)[[1]]
	}

	# Obtain FPKM values for geneList and associated it with gene_ids (outId=T).
	if(doFpkm){
		if(replicates){
			merged <- fpkmRepMat[geneList,]
		}else{
			merged <- fpkmPoolMat[geneList,]
		}
	}else{
		if(replicates){
			merged <- cpmRepMat[geneList,]
		}else{
			merged <- cpmPoolMat[geneList,]
		}
	}	

	# If imposing a mask of samples, use here.
	if(length(mask)>0.1){merged <- .maskSamples(mask,replicates=replicates,fpkmMat=merged)}

	# Populate matrix with gene_id, gene_short_name, and FPKM values.
	merged <- cbind(rownames(merged),idToSym(rownames(merged)),merged)
	colnames(merged)[c(1,2)] <- c('gene_id','gene_short_name')

	# Generate a normalised version of this matrix and extracted melted FPKM values for later.
	normed <- sweep(merged[,-c(1,2)],1,apply(merged[,-c(1,2)],1,max),'/') # Remove 1 and 2 columns (gene info); normalise
	normed <- cbind(merged[,c(1,2)],normed)
	fpkmNorm <- melt(normed,value.var='fpkm',id.vars=c('gene_id','gene_short_name'))$value

	# Melt merged and attached normalised value.
	merged <- melt(merged,value.var='fpkm',id.vars=c('gene_id','gene_short_name'))
	colnames(merged)[c(3,4)] <- c('sample_name','fpkm')
	merged <- cbind(merged,fpkmNorm)

	# Clean up data frame for plotting; change gene_short_name vector to factor in
	# order to manipulate levels for plotting.  The initial rev(.) call plots
	# everything in alphabetical order on the y axis.
	merged$gene_short_name <- as.factor(merged$gene_short_name)

	# Reorder levels for plotting of y axis, in the order received from geneList.
	if(!isGeneId(geneList)){
		merged$gene_short_name <- factor(merged$gene_short_name,levels=geneList)
	}else{
		merged$gene_short_name <- factor(merged$gene_short_name,levels=idToSym(geneList))
	}
		
	preTrunc <- as.data.frame(round(merged$fpkm))
	colnames(preTrunc) <- 'fpkmPre'
	merged <- cbind(merged,preTrunc)

	if (rangeMax<1e5){
		# Track elements that are truncated.
		isTrunc <- merged$fpkm>rangeMax
		merged$fpkm <- pmin(merged$fpkm,rangeMax-0.0001)
	}

	# Define ggplot object.
	heatGr <- ggplot(merged,aes_string(x='gene_short_name',y='sample_name',label='fpkmPre'))

	# Define tile colouring.
	if (doNorm){
		heatGr <- heatGr + geom_tile(aes(fill=fpkmNorm),colour='white')
		theRangeMax <- 1
	}else{
		heatGr <- heatGr + geom_tile(aes(fill=fpkm),colour='white')
		theRangeMax <- max(merged$fpkm)
	}
	if (rangeMax<1e5){
		theRangeMax <- rangeMax+0.0001
	}

	# Define other stylistic features.
	heatGr <- heatGr + scale_fill_gradientn(colours=c('cyan','white','red'),
		values=c(0,0.5,1),guide='colorbar',limits=c(0,theRangeMax)) # Changed from green and purple.
	heatGr <- heatGr + theme(panel.background=element_blank())
	heatGr <- heatGr + theme(axis.text.x = element_text(angle=45,hjust=1))
	heatGr <- heatGr + labs(x='Gene',y='Cell type')

	heatGr <- heatGr + geom_text(data=subset(merged,merged$fpkmPre>rangeMax),aes_string(label='fpkmPre'),size=2)

	# Incorporate any other features.
	heatGr <- heatGr + theme(...)

	print(heatGr)

	invisible(heatGr)
}
