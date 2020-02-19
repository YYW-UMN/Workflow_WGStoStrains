library(ggplot2)
library(plyr)
library(reshape2)
AnalyzeRepeats_dynamic_baseline <- function(farm,max.duplicate,meta.cols) {
	resultlist <- list() 
	for (c in (1:length(unique(farm$CowID)))) {
		cowID <- unique(farm$CowID)[c] # animal ID
		cow <- farm[which(farm$CowID == cowID),] # all data for this animal
                res <- matrix(0,nrow = nrow(cow),ncol = 16)
		cow <- cow[with(cow,order(Days)),] # first, order the data by sampling dates 

                source('/Users/yuanyuanwang/mixed/functions/GetMetaColsIndex.r')
		meta.cols <- GetMetaColsIndex(cow)
		geno <- cow[,-meta.cols]
                snp <- cbind(geno[,which(colSums(geno) > 0 & colSums(geno) < 998)], cow[,meta.cols])

		source('/Users/yuanyuanwang/mixed/functions/CountSNPs.r')
		snp.count <- CountSNPs(snp)

		baseline <- which(snp[1,] > 0)
                baseline.het <- which(snp[1,] >= 10) 

		for (s in (2:nrow(cow))){
			res[s,1] <- cow[s,'State'] 
			res[s,2] <- cow[s,'CowID']
			# total SNPs
			res[s,3] <- rowSums(snp[s,] > 0) 
			# count of homozygous alt (1,2,3,4)
			res[s,4] <- rowSums(snp[s,] >0 & snp[s,] <= 4)
			# count of heterozygous alt (10,20,12,23)
                        res[s,8] <- rowSums(snp[s,] >= 10)
			if (cow[s,'Days'] == cow[s-1,'Days']) { # If sampled on the same day, 
				cow <- cow[with(cow,order('Type2_num')),] 



				
			# count of homozygous alt shared with previous one 
			res[s,5] <- length(intersect(baseline,which(snp[s,] > 0))) # common
	        	res[s,6] <- length(setdiff(which(snp[s,] > 0) ,baseline)) # gained
	        	res[s,7] <- length(setdiff(baseline,which(snp[s,] > 0) )) # lost 
                        res[s,8] <- rowSums(snp[s,] >= 10) # het
			res[s,9] <- length(intersect(baseline.het,which(snp[s,] >= 10))) # common het
			res[s,10] <- length(setdiff(which(snp[s,] >= 10),baseline.het)) # gained het
			res[s,11] <- length(setdiff(baseline.het,which(snp[s,] >= 10))) # lost het
			res[s,12] <- cow[s,'Days'] ; res[s,13] <- cow[s,'Type']; res[s,14] <- cow[s,'Year'];res[s,15] <- cow[s,'Month'];res[s,16] <- cow[s,'Day']
	        }
		resultlist[[c]] <- res
		# plot below
		dat <- cbind(cow[,c('Days','Type')],snp)
		dat$Type <- as.numeric(dat$Type)
		dat$Days <- as.numeric(as.character(dat$Days))
		if (anyDuplicated(dat$Days) > 0) {
			d <- which(duplicated(dat$Days))
			tmp <- seq(1,max.duplicate,1)
			for (i in (1:length(d))) {
				dat$Days[d][i] <- dat$Days[d][i] + tmp[i]
			}
		}
		dat3 <- melt(dat, id.var = c('Days','Type')); colnames(dat3) <- c('Days','Type','SNP','SNP_type')
		dat3$Days = as.factor(dat3$Days)
		snpcolors <- c('0'='cornsilk2','1'='navyblue','2'='red','3'='green','4'='yellow3',
			       '10'='deepskyblue3','20'='salmon','12'='violetred2','23'='burlywood4')	
		if (length(which(dat3$Type == 2)) == nrow(dat3)) {
			na <- rep(1,nrow(dat3)) # all samples were tissue samples 
			p <- ggplot(dat3, aes(SNP, Days)) +
			     geom_tile(aes(fill = as.factor(SNP_type)), colour = 'white') + 
			     scale_fill_manual(values=snpcolors, name = 'SNP type', breaks = c(0,1,2,3,4,10,20,12,23),
					  labels = c('hom ref','hom alt','hom alt2','hom alt3','hom alt4',
						     'het:ref+alt','het:ref+alt2','het:alt1+alt2','het:alt2+alt3')) +
		             theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
			     xlab(paste(cow[1,'State'],' animal ID ',cow[1,'CowID'])) +
			     geom_point(aes(size=ifelse(na,'dot'),color = 'tissue sample')) + 
			     scale_size_manual(values=c(dot=0.45, no_dot=NA), guide="none")
			ggsave(p,file = paste0(cow[1,'State'],'_',cow[1,'CowID'],'.pdf'),scale = 2)
			print(p)
		} else {
			na  <- mapvalues(dat3$Type, c(1,2),c(FALSE,TRUE))
                        p <- ggplot(dat3, aes(SNP, Days)) +
                             geom_tile(aes(fill = as.factor(SNP_type)), colour = 'white') +
                             scale_fill_manual(values=snpcolors, name = 'SNP type', breaks = c(0,1,2,3,4,10,20,12,23),
                                          labels = c('hom ref','hom alt','hom alt2','hom alt3','hom alt4',
                                                     'het:ref+alt','het:ref+alt2','het:alt1+alt2','het:alt2+alt3')) +
                             theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
                             xlab(paste(cow[1,'State'],' animal ID ',cow[1,'CowID'])) +
                             geom_point(aes(size=ifelse(na,'dot','no_dot'),color = 'tissue sample')) +
                             scale_size_manual(values=c(dot=0.45, no_dot=NA), guide="none")
                        ggsave(p,file = paste0(cow[1,'State'],'_',cow[1,'CowID'],'.pdf'),scale = 2)
                        print(p)
 	      	}
	}
	results = do.call(rbind, resultlist)
	colnames(results) <- c('state','cow','total','hom','comm','gained','lost','het','com_het','gai_het','los_het','days','type','year','month','day')
	write.csv(results,file = paste('Results_',farm$State[1],'.csv',sep = ''))
}


























