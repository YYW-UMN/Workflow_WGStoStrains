
RemoveSNPs <- function(gen,counts) {
        g.start <- grep('X',colnames(gen))[1]
        g.end <- tail(grep('X',colnames(gen)), n =1)
	cols <- seq(1,ncol(gen),1) 
	gt.cols <- grep('X',colnames(gen))
	meta.cols <- setdiff(cols,gt.cols)

	dat <- gen[,g.start:g.end]
	remov <- which(counts[12,] == 0)
	dat <- dat[,-remov]
	m <- gen[,meta.cols]
	result <- cbind(m,dat)
	return(result)

}



