GetSNPs <- function(gen){
        snp.pos <- list()
        for (d in 1:length(gen)){
                g <- gen[[d]]
                meta.cols <- GetMetaColsIndex(g)
                str.snp <- colnames(g)[-meta.cols]
                tmp.pos <- rep(0,length(str.snp))
                for (p in (1:length(tmp.pos))){
                        tmp.pos[p] <- as.numeric(strsplit(str.snp[p],'X')[[1]][2])
                }
                snp.pos[[d]] <- tmp.pos
        }
        return(snp.pos)
}

