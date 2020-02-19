SummarizeFarms <- function(farm.data,percent){
        farms <- list()
        farms.geno <- list()
        farms.pos <- list()
        farms.counts <- list()
        farms.summary <- matrix(0, ncol = 9, nrow = length(farm.data))
        colnames(farms.summary) <- c('sequences','animals','repeated','years','totalSNP','homSNP','hetSNP','mixed','rareSNP')
        rownames(farms.summary) <- c(paste('Farm', farm.number,sep = ''))

        for (f in (1:length(farm.data))){
                farm <- read.csv(paste(fpath,farm.data[f],sep = ''),header = TRUE); farm<- farm[,-1] #with both meta data and geno data
                col.start <- which(grepl("X",colnames(farm))==TRUE)[1]-1
                geno <- farm[,-c(1:col.start)] # sample by snp 

                # count SNPs: all SNPs
                pos <- rep(0,length(colnames(geno)))
                for (p in (1:length(pos))){
                        pos[p] <- as.numeric(as.character(substring(colnames(geno)[p],2)))
                }               
                
                counts <- matrix(0,ncol = 4,nrow = ncol(geno))  # rows(snp) by columns(category)
                for (s in (1:ncol(geno))){ 
                        counts[s,1] <- sum(geno[,s] == 1)+ sum(geno[,s] == 2)+ sum(geno[,s] == 3)+ sum(geno[,s] == 4) # hom 
                        counts[s,2] <- sum(geno[,s] == 10)+ sum(geno[,s] == 20)+ sum(geno[,s] == 12)+ sum(geno[,s] == 23) # het
                        counts[s,3] <- sum(geno[,s] == 999) + sum(geno[,s] == 998) # missing value or low coverage
                        counts[s,4] <- sum(geno[,s] != 0) # all SNPs
		}
                
		# remove snp with missing values or low coverage
                rm.snp <- which((counts[,1] == 0) & (counts[,2] == 0))  # meaning that all were missing value or low coverage
		if(length(rm.snp) != 0) {
			 farm <- farm[,-c(rm.snp+col.start)]
			 counts <- counts[-rm.snp,]
		}
                
                farms[[f]] <- farm
                farms.geno[[f]] <- geno
                farms.pos[[f]] <- pos
                farms.counts[[f]] <- counts

                farms.summary[f,1] <- nrow(geno) # total sequences 
                animals <- length(unique(farm[,'CowID'])) # total animals
                farms.summary[f,2] <- animals
                repeats <- which(duplicated(farm[,'CowID']))
                farms.summary[f,3] <- length(which(unique(farm$CowID[repeats])< 7777))# repeated animals
                farms.summary[f,4] <- length(unique(farm[,'Date'])) # years
                farms.summary[f,5] <- nrow(counts) # total SNPs
                farms.summary[f,6] <- sum(counts[,1] != 0 & counts[,2] == 0) # hom SNPs
                farms.summary[f,7] <- sum(counts[,1] == 0 & counts[,2] != 0)  # het SNPs
                farms.summary[f,8] <- sum(counts[,1] != 0 & counts[,2] != 0)  # mix of hom and het SNPs
                farms.summary[f,9] <- sum(counts[,4] < round(animals*percent)) # rare SNPs
	}
	write.csv(farms.summary, file = paste(fpath,'Summarizing_farms.csv',sep = ''))
}
