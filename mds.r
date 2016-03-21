#################################################################################
#
# Mark Cembrowski, Janelia Research Campus, April 3 2015
#
# This is an extension of cummeRbund's MDS that allows for unmasking or masking
# of specific samples. 
#
#################################################################################

# square=F: logical.  Square axes? (ie centred about (0,0)?)
# dims=2: numeric.  Number of dimensions.  Note I've hardcoded this to be either
# 	two or three dimensions.  Any more makes this unwieldly to plot, and 
#	an error is thrown if the user tries to do anything other than these two
#	values.
mds <- function(replicates=F,unmask=c(),mask=c(),square=F,dims=2,x='M1',y='M2',
			doFpkm=T){
	if(doFpkm){
		if(replicates){
			theSamples <- colnames(fpkmRepMat)
			dat <- fpkmRepMat
		}else{
			theSamples <- colnames(fpkmPoolMat)
			dat <- fpkmPoolMat
		}
	}else{
		if(replicates){
			theSamples <- colnames(cpmRepMat)
			dat <- cpmRepMat
		}else{
			theSamples <- colnames(cpmPoolMat)
			dat <- cpmPoolMat
		}
	}

	# Check to make sure mask and unmask are not on at the same time.
	if( ( length(mask)*length(unmask) ) > 0.1){
		stop('Cannot select mask and unmask at same time')
	}

	if(length(mask)>0.1){
		isMasked <- vector(length=length(theSamples))*0
		for (ii in 1:length(mask)){
			isMasked <- isMasked + grepl(mask[ii],theSamples)
		}
		isMasked <- isMasked>0.1
		dat <- dat[,!isMasked]
	}
	if(length(unmask)>0.1){
		isUnmasked <- vector(length=length(theSamples))*0
		for (ii in 1:length(unmask)){
			isUnmasked <- isUnmasked + grepl(unmask[ii],theSamples)
		}
		isUnmasked <- isUnmasked>0.1
		dat <- dat[,isUnmasked]
	}

	dat <- log10(dat+1)

	d <- JSdist(makeprobs(dat))
	fit <- cmdscale(d,eig=F,k=dims)
	if(dims==2){
		res <- data.frame(names=rownames(fit),M1=fit[,1],M2=fit[,2])
		gg <- ggplot(res,aes(color=names,label=names))
	        gg <- gg + geom_point(aes(x=M1,y=M2))
		gg <- gg + geom_text(aes(x=M1,y=M2)) 
	}else{
		if(dims!=3){
			stop('Need either two or three dimensions.')
		}else{
			res <- data.frame(names=rownames(fit),M1=fit[,1],M2=fit[,2],M3=fit[,3])
			gg <- ggplot(res,aes(color=names,label=names))
			gg <- gg + geom_point(aes_string(x=x,y=y)) 
			gg <- gg + geom_text(aes_string(x=x,y=y))
		}
	}

	maxAxis <- max(fit)*1.1
	gg <- gg + expand_limits(x=c(-maxAxis,maxAxis),y=c(-maxAxis,maxAxis))
	gg <- gg + theme_bw()
	gg <- gg + theme(legend.position='none')
	print(gg)


#	if(square){
#		temp <- ggplot_build(gg)
#		xMax <- max(abs(temp$panel$ranges[[1]]$x.range))
#		yMax <- max(abs(temp$panel$ranges[[1]]$y.range))
#		gg <- gg + xlim(-xMax,xMax)
#		gg <- gg + ylim(-yMax,yMax)
#	}
}
	
