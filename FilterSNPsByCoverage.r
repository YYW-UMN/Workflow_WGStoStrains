library(dplyr)
library(plyr)

FilterSNPsByCoverage <- function(low.cov,gt,dp) {
# 4. Filtering out SNPs with low coverage across samples
  low.cov <- 10
  gtt <- gt %>% mutate_all(as.character)
  gtt <- as.matrix(gtt)
  gtt[which(dp < low.cov)] <- '999'
  GT <- mapvalues(gtt, c('0/0','1/1','2/2','3/3','4/4','0/1','0/2','1/2','2/3','.','999'), c(0,1,2,3,4,10,20,12,23,999,998))
  class(GT) <- "numeric"
  write.csv(GT,paste(fpath,'GT_LowCoverageLabeled.csv', sep = ''))
}



