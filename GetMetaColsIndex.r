GetMetaColsIndex <- function(gen){
        cols <- seq(1,ncol(gen),1)
        gt.cols <- grep('X',colnames(gen))
        meta.cols <- setdiff(cols,gt.cols)
        return(meta.cols)
}

