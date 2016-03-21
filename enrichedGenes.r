###################################################################################################
#
# Mark Cembrowski, Janelia Farm, October 1 2015
#
# This script looks for enriched (or conversely, depleted) genes in a given
# sample or samples, relative to all other samples in the dataset.  
#
# Note that this script requires ALL replicates to be at least X-fold enriched
# relative to ALL OTHER replicates.  This is stringent, but helps to eliminate
# outliers.  This can be circumvented by using the avgPass=T, which operates
# on the average values only.
#
# 111213: EDIT: it would be advantageous to be able to calculate enriched genes while masking
# some replicates from the comparison.  I've implemented a mask option to do this.
#
# 100115: EDIT: I've added functionality that allows enriched genes to be obtained from CPM,
# rather than FPKM, data.
#
# INPUT:
# sampleNames: a vector of one or more sample names.
# OPTIONAL ARGUMENTS:
# foldThres: fold difference to count as enriched [default=3]
# fpkmThres: require a threshold for minimum FPKM value in enriched popln [default=-1; ie off]
# depleted: search for depleted genes within this population instead? [default=F]
# avgPass: relax pass criteria by considering only averages (see header of this
# 	function for explanation)?  [default=F]
# mask: a vector of one or more sample names to mask from analysis (ie, they will not act in
#	determining the enrichment cutoff).  [default=c(); ie off]
# plotAvg: plot a graph of the average fold change between samples?  [default=F]
#
# OUTPUT:
# A list is returned with all XLOC value that fulfill the given criteria.
#
####################################################################################################

enrichedGenes <- function(sampleNames,foldThres=3,fpkmThres=-1,depleted=F,avgPass=F,mask=c(),plotAvg=F,shuffled=F,doFpkm=T){
	# Generate FPKM matrix to operate on.
	if(doFpkm){
		if(avgPass){
			theMat <- fpkmPoolMat
		}else{
			theMat <- fpkmRepMat
		}
	}else{
		if(avgPass){
			theMat <- cpmPoolMat
		}else{
			theMat <- cpmRepMat
		}
	}
		

	toMask <- c() # Set up a dummy variable that will remain blank if not using mask.
	
	# Identify replicates to mask, if selected.
		if(length(mask)>0.1){
			toMask <- grepl(paste(mask,collapse='|'),colnames(theMat))
			theMat <- theMat[,!toMask]
		}

	if (avgPass){
		sampleReps <- colnames(theMat)%in%sampleNames
	}else{
		# Extract features of samples and replicates.
		reps <- colnames(theMat)
		samples <- unique(substr(reps,1,nchar(reps)-2)) # Remove _X suffix from reps
		
		# Verify that supplied sample names are all from dataset.
		if( prod(sampleNames%in%samples)<0.1 ){
			stop('At least one supplied sample not in data set sample names.')
		}
			
		# Identify replicates that are part of samples to identify fold changes for.
		sampleReps <- ( (substr(reps,1,nchar(reps)-2) )%in%sampleNames)

	}
	
	if(shuffled){
		theMat <- .shuffleMat(theMat)
	}
	
	
	# Check if doing enrichment or depletion survey; proceed accordingly.
	if(!depleted){
		# Get minimum replicate FPKM in supplied population.
		minFpkm <- apply(as.matrix(theMat[,sampleReps]),1,min) #as.matrix prevents an error when theMat[,sampleReps] has one column and becomes a vector

		# Get maximum replicate FPKM in rest of population.
		maxFpkm <- apply(as.matrix(theMat[,!sampleReps]),1,max)
	
		# Mask FPKM values below cutoff.
		maskFpkm <- minFpkm>fpkmThres
		minFpkm <- minFpkm*maskFpkm
		
	}else{
		# Get maximum replicate FPKM in supplied population.
		maxFpkm <- apply(as.matrix(theMat[,sampleReps]),1,max)
		# Get minimum replicate FPM in rest of population.
		minFpkm <- apply(as.matrix(theMat[,!sampleReps]),1,min)
		
		# Mask FPKM values below cutoff.
		maskFpkm <- minFpkm>fpkmThres
		minFpkm <- minFpkm*maskFpkm
	}

	# Check if fold difference surpassed.
	theOut <- rownames(theMat)[ minFpkm > (maxFpkm*foldThres) ]

	# Plot graph illustrating fold change, if desired.
	if(plotAvg){
		# Reduce matrix to those entries crossing fold difference.
		plotMat <- theMat[ minFpkm > (maxFpkm*foldThres) , ]
		
		# From the above code, columns have already been assigned to being desired sample
		# reps (sampleReps) and masked samples (toMask).  Use this to generate the
		# mean FPKM values for desired samples and comparison samples, masking any
		# contributions from undesired comparisons.  Note that these have slightly different
		# formats -- sampleReps is a string of logicals equal to the numbe of columns, 
		# whereas toMask gives the indices of masked columns.  This requires a 
		# housekeeping step.
		#
		# In the below, I'll refer to the supplied samples and comparison samples as
		# pop1 and pop2, respectively -- this general naming is chosen because they could
		# corresponding to either the enriched or depleted sample, set by the 'depleted'
		# option.
		pop1 <- apply(as.matrix(plotMat[,sampleReps]),1,mean)
		sampleRepsToInd <- which(sampleReps>0.9)
		pop2 <- apply(as.matrix(plotMat[,-sort(c(sampleRepsToInd,toMask))]),1,mean)
		popEnr <- pmax(pop1/pop2,pop2/pop1) # Fold change, take >1
		compPops <- as.data.frame(cbind(pop1,pop2,popEnr))
		
		thePlot <- ggplot(compPops,aes(factor(pop1*0),popEnr)) 
		thePlot <- thePlot + geom_boxplot(outlier.size=0) + geom_jitter() # geom_jitter includes all outliers;
							 			  # otherwise left with two data points
										  # for each outlier.
		thePlot <- thePlot + expand_limits(y=0)
		
		print(thePlot)
		
		#return(compPops)
	}

	invisible(theOut)
}

# Shuffle a matrix.  Used for null hypothesis testing.
.shuffleMat <- function(theMat){
	trans <- t(apply(theMat,1,sample))
	colnames(trans) <- colnames(theMat)
	invisible(trans)
} 
