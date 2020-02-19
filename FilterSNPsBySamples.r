
FilterSNPsBySamples <- function(percentage,cutoff, gt_labeled,dp) {
        snps <- gt_labeled
        remov <- c()
        cf <- cutoff
        counts <- matrix(0,nrow = nrow(snps),ncol = 2) # two columns: nonzero and 999/998
        for (i in (1:nrow(snps))) {
                counts[i,1] <- sum(snps[i,] != 0)
                counts[i,2] <- sum(snps[i,] == 999) + sum(snps[i,] == 998) 
                if ((counts[i,1] - counts[i,2]) < cf) {
                        remov <- append(remov,i)
                }
        }
        snps_f <- snps[-remov,]
        dp_f <- dp[-remov,]
        cat('Using a cutoff at',cf,',',length(remov),'rare SNPs were removed from a total of',nrow(snps),'\n', sep = ' ') 

	write.csv(snps_f,file = paste(fpath,'GT_Filtered_AtLeast_',cf,'_Animals.csv', sep = ''))
        write.csv(dp_f,file = paste(fpath,'DP_Filtered_AtLeast_',cf,'_Animals.csv', sep = ''))
        write.csv(remov,file = paste(fpath,'Removed_Positions_AtLeast_',cf,'_Animals.csv', sep = ''))
}

