CountSNPsByVar <- function(dat.list) {
        count.list <- list() # store counts 
        updat.list <- list() # store updated data (zero SNPs removed)
        final.count.list <- list() # count again 
        for (c in (1:length(dat.list))){
                count.list[[c]] <- CountSNPs(dat.list[[c]])
                updat.list[[c]] <- RemoveSNPs(dat.list[[c]],count.list[[c]])
                final.count.list[[c]] <- CountSNPs(updat.list[[c]])
        }
        return(final.count.list)
} 

