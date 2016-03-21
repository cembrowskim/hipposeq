###################################################################
#
# Mark Cembrowski, Janelia Farm, Jan 29 2014
#
# This function takes a string of gene symbols and returns all
# associated gene IDs.
#
# INPUT:
# syms: character.  one or more gene symbols to convert.
#
# OPTIONAL ARGUMENTS:
# keepMissed=F: logical.  if supplied with symbol(s) that do not
# match any gene_short_name, keep these as 'NA' in the output?
# quiet=F: logical. suppress printout of any missed conversions.
#
###################################################################

symToId <- function(syms,keepMissed=F,quiet=F){

	# Load up vectors corresponding to gene short names and
	# gene ids.
	numSamps <- length(unique(fpkmPool$sample_name))
	numInds <- nrow(fpkmPool)/numSamps
	allGeneIds <- as.vector( (fpkmPool$gene_id)[1:numInds] )
	allGsns <- as.vector( (fpkmPool$gene_short_name)[1:numInds] )

	# Identify indices associated with supplied short names, and
	# convert to gene IDs.
	theInds <- which(allGsns%in%syms)
	theIds <- allGeneIds[theInds]
	names(theIds) <- allGsns[theInds]

	theIds <- theIds[syms] # Order according to order of input symbols.

	# Screen for any missed symbols; announce.
	missed <- syms[!(syms%in%allGsns)]

	if(length(missed)>0.1){
		if(!quiet){
			print("Note: missed the following during a sym -> id conversion:")
			print(missed)
		}
	
		# Remove from list, if option selected.
		if(!keepMissed){
			theIds <- theIds[-which(!(syms%in%allGsns))]
		}
	}

	theOut <- list(theIds,missed)

	invisible(theOut)
}
