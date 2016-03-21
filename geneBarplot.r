#####################################################################################
#
# Mark Cembrowski, Janelia Farm, January 28 2014
#
# The following script returns a barplot of expression for a given gene, with
# confidence intervals as error bars.
#
# 1011515: updated to include options on CPM-based, rather than FPKM-based, data.
#
# INPUTS:
#  theGene: character.  Name of gene to plot.
# OPTIONAL INPUTS:
#  replicates=F: logical.  Display data points of replicates?
#  mask=c(): character.  Mask particular samples from being plotted?  Should not be
#	used in conjunction with unmask option.
#  unmask=c(): character.  Include only particular samples on plot?  Should not be
# 	used in conjunction with mask option.
#  rMax=-1: variable.  Send the max range on the graph.
#
#####################################################################################

geneBarplot <- function(theGene,replicates=F,mask=c(),unmask=c(),rMax=-1,doFpkm=T){
	if(length(mask)*length(unmask)>0.1){
		stop('Cannot use mask and unmask options at same time.')
	}

	if(!isGeneId(theGene)){
		# Transform to Ensembl ID.
		theGene <- symToId(theGene)[[1]] # Successful hits.
	}

	# Get FPKM value for gene.
	if(doFpkm){
		df <- subset(fpkmPool,gene_id==theGene)
	}else{
		df <- subset(cpmPool,gene_id==theGene)
	}

	# Get replicates.
	if(doFpkm){
		dfRep <- subset(fpkmRep,gene_id==theGene)
	}else{
		dfRep <- subset(cpmRep,gene_id==theGene)
	}

	# Mask samples, if selected.
	if(length(mask)>0.1){
		toMask <- grep(paste(mask,collapse='|'),df$sample_name)
		if(length(toMask)<0.1){
			stop('All requested masked sample(s) not found.')
		}
		
		df <- df[-toMask,]

		if(replicates){
			toMask <- grep(paste(mask,collapse='|'),dfRep$sample_name)
			dfRep <- dfRep[-toMask,]
		}
	}	

	# Show only unmasked samples, if selected.  This operation is the opposite
	# of mask; ie, only the unmasked samples will be shown.
	if(length(unmask)>0.1){
		toUnmask <- grep(paste(unmask,collapse='|'),df$sample_name)
		if(length(toUnmask)<0.1){
			stop('All requested unmasked sample(s) not found.')
		}

		df <- df[toUnmask,]

		if(replicates){
			toUnmask <- grep(paste(unmask,collapse='|'),dfRep$sample_name)
			dfRep <- dfRep[toUnmask,]
		}
	}

	# Generalise column 'fpkm' or 'cpm' name to 'metric'.
	if(doFpkm){
		metricColPool <- grep('fpkm',colnames(df))
		metricColRep <- grep('fpkm',colnames(dfRep))
	}else{
		metricColPool <- grep('cpm',colnames(df))
		metricColRep <- grep('cpm',colnames(dfRep))
	}
	colnames(df)[metricColPool] <- 'metric'
	colnames(dfRep)[metricColRep] <- 'metric'

	# Plot.
	p <- ggplot(df,aes(x=sample_name,y=metric,fill=sample_name))
	p <- p + geom_bar(stat='identity')
	if(replicates){
		p <- p + geom_point(aes(x=sample_name,y=metric),size=3,
			shape=18,colour='black',data=dfRep)
	}

	if(doFpkm){
		p <- p + geom_errorbar(aes(ymin=conf_lo,ymax=conf_hi,group=1),width=0.5)
	}else{
		# No CI for CPM.
	}

	if(rMax>0){
		p <- p + ylim(c(0,rMax))
	}

	p <- p + xlab('Cell type') 
	if(doFpkm){
		p <- p + ylab('FPKM')
	}else{
		p <- p + ylab('CPM')
	}
	p <- p + ggtitle(idToSym(theGene))
	p <- p + theme_bw() + theme(legend.position='none')

	print(p)
}	
