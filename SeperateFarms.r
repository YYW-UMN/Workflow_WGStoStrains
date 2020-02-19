library(plyr)
library(dplyr)


  SeperateFarms <- function(meta_data,genotype_data) {
        ma <- meta_data # dim(samples info)
        ga <- genotype_data # dim(snps sample)

        fm <- unique(ma$Farm)
        fnames <- do.call("paste", c(list(rep('Farm',length(fm)),fm), sep = ''))
        farms <- list()

        for (f in (1:length(fm))){
                farm.id <- unique(ma$Farm)[f]
                samples <- which(ma$Farm == farm.id) # samples
                snps <- which(rowSums(ga[,frows])>0) # non-zero snps 
                ma.tmp <- ma[samples,]
                ga.tmp <- ga[snps,samples]
                farms[[f]] <- cbind(ma.tmp,t(ga.tmp))
                write.csv(file = paste(fpath,fnames[f],'_SampleBySNP.csv',sep = ''), farms[[f]])
        }
  }


