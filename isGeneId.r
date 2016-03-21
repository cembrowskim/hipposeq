##############################################################################
#
# Mark Cembrowski, Janelia Farm, October 26 2013
#
# This script returns a logical, based upon whether the element/list
# contains gene_ids.  This is a quick way of differentiating between lists
# of gene_short_names and lists of gene_ids.  Note that gene_ids can take
# the form of XLOC values (if transcriptome reconstructed de novo) or
# annotated transcriptome values (e.g., ENSMUSG...) if mapped using a 
# reference transcriptome.
#
# This is a convenience function used by many scripts. Input/output self
# explanatory.
#
#############################################################################

isGeneId <- function(geneList){
	if(grepl('XLOC|ENSMUS',geneList[1])){
		aGeneId <- T
	}else{
		aGeneId <- F
	}
	invisible(aGeneId)
}
