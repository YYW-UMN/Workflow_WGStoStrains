UpdateData <- function(dat.list) {
        count.list <- list() # store counts 
        updat.list <- list() # store updated data (zero SNPs removed)
        for (c in (1:length(dat.list))){
                count.list[[c]] <- CountSNPs(dat.list[[c]])
                updat.list[[c]] <- RemoveSNPs(dat.list[[c]],count.list[[c]])
        }
        return(updat.list)
}

