#####################################################################################################
#
# Mark Cembrowski, Janelia Research Campus, April 25 2015
#
# This script plots a scatter of count values.
#
# It is analogous to the fpkmScatter approach used for FPKM data.
#
# OPTIONAL ARGUMENTS:
# cpm=T: logical.  Use counts per million (default) or raw counts?
# quiet=T: logical. If F, print out correlation coefficient and number of DE genes.
#
#####################################################################################################

geneScatter <- function(sampleX,sampleY,fdr=0.05,doFpkm=T,maxX=-1,minX=-1,minY=-1,maxY=-1,cpm=T,
			label=F,colour=c('green'),gsnList=c(),quiet=T){

	# Check if supplied a replicate, which means suffix will be _X, where
	# X is an integer.
	suff1 <- substr(sampleX,nchar(sampleX)-1,nchar(sampleX)-1) 
	suff2 <- substr(sampleX,nchar(sampleX),nchar(sampleX))

	if( (suff1=="_") && (!is.na(suppressWarnings(as.numeric(suff2)))) ){
		doRep <- T
		print('Detected a replicate sample.')
	}else{
		doRep <- F
	}

	# Retrieve values for visualisation and analysis.
	if(doFpkm){
		ids <- rownames(fpkmPoolMat) #same for rep or pool
		if(doRep){
			xVals <- fpkmRepMat[,sampleX]+1 
			yVals <- fpkmRepMat[,sampleY]+1
		}else{
			xVals <- fpkmPoolMat[,sampleX]+1 # additive smoothing
			yVals <- fpkmPoolMat[,sampleY]+1 # additive smoothing
		}
	}else{	
		ids <- rownames(cpmPoolMat) # same for rep or pool; cpm or count
		if(cpm){
			if(doRep){
				xVals <- cpmRepMat[,sampleX]+1
				yVals <- cpmRepMat[,sampleY]+1
			}else{
				xVals <- cpmPoolMat[,sampleX]+1 # additive smoothing
				yVals <- cpmPoolMat[,sampleY]+1 # additive smoothing
			}
		}else{
			if(doRep){
				xVals <- countRepMat[,sampleX]+1
				yVals <- countRepMat[,sampleY]+1
			}else{
				xVals <- countPoolMat[,sampleX]+1 # additive smoothing
				yVals <- countPoolMat[,sampleY]+1 # additive smoothing
			}
		}
	}
	
	# Get DE'ed genes.
	if(length(gsnList)<0.1){
		if(doFpkm){
			isDe <- vector(length=nrow(fpkmPoolMat))>0.1
			if(!doRep){
				names(isDe) <- rownames(fpkmPoolMat)
				theDe <- getSig(cuff_data,sampleX,sampleY,alpha=fdr,level='genes')
				isDe[theDe] <- T
			}
		}else{
			if(!doRep){
				toKeep <- grepl(paste(c(sampleX,sampleY),collapse='|'),colnames(countRepMat))
				subMat <- countRepMat[,toKeep]
		
				toKeep <- grepl(paste(c(sampleX,sampleY),collapse='|'),countRepData$condition)
				subCol <- countRepData[toKeep,]
				subCol$condition <- droplevels(subCol$condition)

				dds <- DESeqDataSetFromMatrix(countData=subMat,colData=subCol,design=~condition)
	
				dds <- DESeq(dds)
				res <- results(dds)

				isDe <- res[,'padj']<fdr
				isDe[is.na(isDe)] <- F
			}else{
				isDe <- vector(length=length(ids))>0.1 # No DE genes for replicates.
			}
		}
	}else{
		if(doFpkm){
			isDe <- vector(length=nrow(fpkmPoolMat))>0.1
			names(isDe) <- rownames(fpkmPoolMat)
		}else{
			isDe <- vector(length=nrow(cpmPoolMat))>0.1
			names(isDe) <- rownames(cpmPoolMat)
		}
		
		if(!isGeneId(gsnList)){
			isDe[unlist(symToId(gsnList))] <- T
		}else{
			isDe[gsnList] <- T
		}
	}


	df <- data.frame(gene_id=ids,valX=xVals,valY=yVals)
	df[,'gene_short_name'] <- unlist(idToSym(ids))
	rownames(df) <- df$gene_id
	
	if(!quiet){
		print(paste("Number of DE genes:",sum(isDe)))
		print(paste("Pearson corr of non-log-transformed vals:",cor(df$valX,df$valY)))
	}

	if(length(colour)>1.1){
		# Have suggested two colours; add in two levels of highlighting values.
		yEnr <- df$valY>df$valX
		yEnrDe <- (yEnr*isDe)>0.1
		toHl <- isDe
		toHl[yEnrDe] <- 2
		df[,'highlight'] <- toHl
	}else{
		df[,'highlight'] <- isDe
	}
	theMax <- max(max(df[,2]),max(df[,3]))
	theMin <- 1
	theMaxX <- theMax
	theMaxY <- theMax
	theMinX <- theMin
	theMinY <- theMin
	if (maxX>0){theMaxX <- maxX}
	if (maxY>0){theMaxY <- maxY}
	if (minX>0){if (minX<maxX){theMinX <- minX}}
	if (minY>0){if (minY<maxY){theMinY <- minY}}
	
	gg <- ggplot(df,aes(x=valX,y=valY,label=gene_short_name)) 
	gg <- gg + geom_point(aes(colour=factor(highlight)),alpha=0.2)
	gg <- gg + scale_colour_manual(values=c("black",colour))
	gg <- gg + theme_bw()
	gg <- gg + scale_x_log10(limits=c(theMin,theMax))
	gg <- gg + scale_y_log10(limits=c(theMin,theMax))
	if(doFpkm){
			gg <- gg + xlab(paste('FPKM+1',sampleX))
			gg <- gg + ylab(paste('FPKM+1',sampleY))
			gg <- gg + geom_abline(slope=1,linetype='dashed',weight=0.1)
	}else{
		if(cpm){
			gg <- gg + xlab(paste('CPM+1,',sampleX))
			gg <- gg + ylab(paste('CPM+1,',sampleY))
			gg <- gg + geom_abline(slope=1,linetype='dashed',weight=0.1)
		}else{
			gg <- gg + xlab(paste('Count+1,',sampleX))
			gg <- gg + ylab(paste('Count+1,',sampleY))
		}
	}
	gg <- gg + theme(legend.position='none')
	if(label){gg <- gg + geom_text(data=subset(df,highlight>0.1),hjust=0,vjust=0,size=2)}


	print(gg)
}

