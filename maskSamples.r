############################################################################
#
# Mark Cembrowski, Janelia Farm, Jan 29 2014
#
# This function is a helper function and returns the FPKM matrix with the
# desired samples masked.
#
# fpkmMat can be supplied if the full fpkmMatrix (or repFpkmMatrix) should
# not be used or if regeneration is undesired (say, for large matrices).
#
# The masking is implemented with grep, meaning it is INCLUSIVE OF ALL
# STRINGS that exhibit partial matching to the mask.
#
# INPUT:
# providedNames: character.  Name of sample(s) to mask.  Masked via grep,
#	not equality.
#
# OPTIONAL ARGUMENTS:
# replicates=F: logical. use replicate data.
# fpkmMat=matrix(): matrix. Can provided specific matrixd for masking; otherwise
#	defaults to full matrix.
#
############################################################################

.maskSamples <- function(providedNames,replicates=F,fpkmMat=matrix()){

	if(length(providedNames)>0.1){
		# Screen to make sure masked elements are present in the
		# sample names, and if so, mask accordingly.
		# Check if samples are prevent.  theMat is the invisibly returned
		# FPKM matrix.
		theMat <- .screenSamples(providedNames,replicates,fpkmMat)

		# Identify column indices to mask.
		strMask <- paste(providedNames,collapse='|')
		colMask <- grep(strMask,names(theMat))

		# Mask.
		theMat <- theMat[,-colMask]

	}else{
		# Just return the matrix that would be used for the
		# ensuing analysis.
		warning('No mask provided to mask call.  Echoing input as output.')
		if(prod(dim(fpkmMat))>1.1){
			theMat <- fpkmMat
		}else{
			if(replicates){
				theMat <- repFpkmMatrix(genes(cuff_data))
			}else{
				theMat <- fpkmMatrix(genes(cuff_data))
			}
		}
	}
	
	invisible(theMat)
}	
