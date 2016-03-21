#################################################################################
#
# Mark Cembrowski, Janelia Farm, January 29 2014
#
# This script converts a list of gene IDs to gene symbols. 
#
# INPUT:
# geneIds: character.  One or more gene IDs to convert.
#
# OUTPUT:
# theGsns: character.  List of symbols corresponding to gene IDs.
#
# TODO: incorporate an alert in case one or more gene IDs are missed during 
# conversion.
#
#################################################################################

idToSym <- function(ids){
	# Load up vectors correspond to gene short names and gene ids.
	numSamps <- length(unique(fpkmPool$sample_name))
	numInds <- nrow(fpkmPool)/numSamps
	allGeneIds <- as.vector( (fpkmPool$gene_id)[1:numInds] )
	allGsns <- as.vector( (fpkmPool$gene_short_name)[1:numInds] )

	# Identify indices associated with supplied IDs, and convert to short
	# names.
	theInds <- which(allGeneIds%in%ids)
	theGsns <- allGsns[theInds]
	names(theGsns) <- allGeneIds[theInds]
	theGsns <- theGsns[ids]
	
	invisible(theGsns)
}
