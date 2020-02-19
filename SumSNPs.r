SumSNPs <- function(gen) {
        m.cols <- GetMetaColsIndex(gen)
        dat.mat <- gen[,-m.cols]
        total.snp <- rowSums((dat.mat != 0)&(dat.mat < 998))
        het.snp <- rowSums((dat.mat > 9)&(dat.mat < 998))
        hom.snp <- rowSums((dat.mat != 0)&(dat.mat < 9))
        sum.totals <- list(total.snp,het.snp,hom.snp)
        return(sum.totals)
}


