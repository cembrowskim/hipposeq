###############################################################################
#
# Mark Cembrowski, Janelia Farm, October 16 2013
#
# UPDATE: this now includes a family of functions for saving ENSEMBL IDs for
# commonly-used lists of genes; ie, DE'ed genes and marker genes.
#
# This function converts a vector of either gene_ids or gene_short_names
# into the corresponding Ensembl IDs, suitable for use in DAVID.
#
# INPUT: a vector of either gene_ids or gene_short_names
# OPTIONAL INPUT:
#	saveName: add the file name if want to write this to a table.  The
#		.txt file extension is added automatically.  This is saved
#		in the current working directory.  [default=''; ie no save]
#	saveDir: add the directory to where this is saved [default=''; ie
#		saved in current working directory].  Currently this is
#		implemented to save in subfolders of the current working
#		directory only.  saveDir is the subdirectory to save to; if
#		the directory exists, the file is added; if the directory
#		does not exist it is created and saved there.
#	parseComma: with cuffmerge, loci that overlap are merged into the 
#		same transcript.  This produces gene names like
#		'Mir5129,Zeb2', which can be deparsed and retrieved prior
#		to saving the gene names.  Do this?  [default=T]
# OUTPUT: a vector of Ensembl IDs.  If a saveName is specified, a text file
#		is also written.
#
###############################################################################

getEnsId <- function(toConvert,saveName='',saveDir='',parseComma=T){
	# First, check if have either gene_ids, gene_short_names, or Ensembl
	# identifiers.  If have gene_ids, convert to short names, and then
	# to Ensembl IDs.  If have short names, convert to Ensembl IDs.  If
	# have Ensembl IDs, then already good to do.
	if(length(grep('XLOC',toConvert[1]))>0.1){
		# Have XLOC values; ie, gene_ids.  Convert to gene_short_names
		print('Interpreting gene list as XLOC values.')
		toConvert <- xlocToSym(toConvert)
	}

	# Execute this block of code if do not already have Ensembl IDs, in
	# order to obtain them.
	if(!(length(grep('ENSMUS',toConvert[1]))>0.1)){
		# Next, convert gene_short_names to Ensembl IDs.  Need to load
		# a lookup table.  For now, this is just comprised of loading a set
		# of Cufflinks data built on the Ensembl transcriptome.
		symbolAndEns <- read.table('~/research/seq/seqData/ca1/ca1_12/genes.fpkm_tracking',header=T,sep='\t')
		symbolAndEns <- symbolAndEns[,c(1,5)] # tracking_id and g_s_n
	
		ensInds <- which(symbolAndEns$gene_short_name%in%toConvert)
		converted <- as.vector(symbolAndEns$tracking_id[ensInds])
	
		# Declare which genes were missed.
		print('The following are nonhits:')
		notIn <- toConvert[!(toConvert%in%symbolAndEns$gene_short_name)]
		print(notIn)
	
		# If deparsing option selected, deparse.
		if(parseComma){
			deparsed <- unlist(strsplit(notIn,','))
			# Remove common sources of merged loci that are typically not
			# on: Mir, Gm, and SNORD genes.
			if(length(grep('Mir|SNORD|Gm',deparsed))>0.1){
				deparsed <- deparsed[-(grep('Mir|SNORD|Gm',deparsed))]
			}
	
			# Let user know what has been saved.
			print('The following were recovered from comma deparsing:')
			print(deparsed)
			
			# Retrieve gene names	
			ensIdsD <- which(symbolAndEns$gene_short_name%in%deparsed)
			convertedD <- symbolAndEns$tracking_id[ensIdsD]
	
			# Add to converted.
			converted <- c(converted,as.vector(convertedD))
		}
	}else{
		converted <- toConvert # toConvert is already Ensembl IDs
	}
	
	# Save the file in the current working directory, if specified
	if(nchar(saveName)>0.5){
		if(nchar(saveDir)>0.5){
			# Create directory if does not exist; ignore warning if it does.
			dir.create(file.path(getwd(),saveDir),showWarnings=F)
			write.table(converted,quote=F,row.names=F,col.names=F,file=paste(saveDir,'/',saveName,'.txt',sep=''))
		}else{
			write.table(converted,quote=F,row.names=F,col.names=F,file=paste(saveName,'.txt',sep=''))
		}
	}

	invisible(as.vector(converted))
}

# The below function extracts ENSEMBL IDs for genes that are commonly needed,
# ie, DE'ed genes and marker genes.  Below, criteria is the argumnent that
# tells the function what criteria to screen for.  criteria='de' is DE'ed
# genes, criteria='mark' is marker genes.  If not entered, the default
# is to return ENSEMBL IDs for DE'ed genes.
#
# This, by default, saves the files in the current working directory, with 
# informative names.  A run of this script will produce three outputs:
# 1.  A list of the ENSEMBL tags for ALL 'criteria' genes are returned [ie,
#	indiscriminately in which population they are enriched in]
# 2 and 3.  Two lists, separated according to the population they are 
# 	enriched in.
#
# OPTIONAL ARGUMENTS:
# theFoldThres and theFpkmThres are paraeters used for extracting marker genes.
standardEnsId <- function(sampleX,sampleY,criteria='',theFoldThres=5,theFpkmThres=5){
	
	defaultCriteria <- 0
	if(length(grep('mark',criteria))<0.1){
		if(length(grep('de',criteria))<0.1){
			print('No valid criteria entered for selecting a list of genes to return ENSEMBL IDs.  Defaulting to DE genes.')
			defaultCriteria <- 1
		}
	}

	# Generate geneList for the proper criteria and assign save name.
	if( (length(grep('de',criteria))>0.9) || (defaultCriteria>0.9) ){
		# Obtain DE'ed genes.
		genToUse <- getSigTable(cuff_data,alpha=0.05,level='genes')

		# If have >2 samples, make sure both are in sigTable.  [If
		# have only two samples, getSigTable returns a vector; if
		# >2, returns a matrix.]  Once this check is completed,
		# write the gene list of interest.
		compName <- paste(sampleX,sampleY,sep='vs')
		if(!is.vector(genToUse)){
			# Check to make sure both samples are in sigTable
			indOfComparison <-0 # index of column containing relevant comparison
			indOfComparison <- grep(compName,colnames(genToUse)) # index of column containing relevant comparison
			if(length(indOfComparison)<0.1){
				# The proper index will not be returned if the sample
				# names are juxtaposed.  Try opposite order.
				compName <- paste(sampleY,sampleX,sep='vs')
				indOfComparison <- grep(compName,colnames(genToUse))
				if(length(indOfComparison)<0.1){
					print('Could not locate relevant sample comparison.')
					print('Check to make sure both entered correctly.')
					print('Aborting.')
					return
				}
			}
	

			# Remove NA values and zero values; map to XLOC values.
			geneList <- genToUse[,indOfComparison]
			geneList[is.na(geneList)] <- 0
			geneList <- names(geneList[geneList>0.1])
			
		}else{
			geneList <- names(genToUse)
		}

		# Next, need to break this down into lists depending on
		# whether the gene is enriched in sample X or Y.
		detailed <- pairwiseDetails(sampleX,sampleY,geneList)
		enrX <- rownames(detailed)[detailed$overexpressedIn=='x']
		enrY <- rownames(detailed)[detailed$overexpressedIn=='y']
		

		saveNameAll <- paste('ensDeGenes',compName,sep='')
		saveNameXEnr <- paste('ensDeGenes',sampleX,'EnrOver',sampleY,sep='')
		saveNameYEnr <- paste('ensDeGenes',sampleY,'EnrOver',sampleX,sep='')
	}else{
		if(length(grep('mark',criteria))>0.9){
			geneList <- foldDiffThres(sampleX,sampleY,foldThres=theFoldThres,fpkmThres=theFpkmThres)
			enrX <- as.vector(geneList$gene_short_name[geneList$overexpressedIn=='x'])
			enrY <- as.vector(geneList$gene_short_name[geneList$overexpressedIn=='y'])
			geneList <- as.vector(geneList$gene_short_name)

			saveNameAll <- paste('ensMarkerGenes',sampleX,'vs',sampleY,sep='')
			saveNameXEnr <- paste('ensMarkerGenes',sampleX,'EnrOver',sampleY,sep='')
			saveNameYEnr <- paste('ensMarkerGenes',sampleY,'EnrOver',sampleX,sep='')
		}else{
			print('hm.  seemed to have entered an invalid criteria.')
			print('aborting. step up your game')
			return()
		}
	}
	print('Doing full list...')
	getEnsId(geneList,saveName=saveNameAll)
	print('Doing x-enriched list...')
	getEnsId(enrX,saveName=saveNameXEnr)
	print('Doing y-enriched list...')
	getEnsId(enrY,saveName=saveNameYEnr)
}
			
