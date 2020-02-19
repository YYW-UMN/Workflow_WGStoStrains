library(apcluster)
library(reshape)
library(scales)
source('GetMetaColsIndex.r')

gt <- read.csv('~/mixed/output/GT_2755_sample_by_snp_with_meta_cols.csv', header = TRUE); gt <- gt[,-c(1,2)]
meta.cols <- GetMetaColsIndex(gt)
meta <- gt[,meta.cols]
meta$idx <- seq(1,nrow(meta),1)
g <- gt[,-meta.cols]
ann <- which(meta$Contamination == 1)

# similarity matrix
sim <- negDistMat(g)
msim <- melt(sim)

# first add the farm label for the opponent in the comparison
colnames(msim) <- c('idx','itself','similarity')
msim <- join(msim,meta[,c('idx','Farm')],by = 'idx')

# second add the farm label for the itself in the comparison
colnames(msim) <- c('opponent','idx','similarity','oppo.farm')
msim <- join(msim,meta[,c('idx','Farm')],by = 'idx')

# finalize the column names
colnames(msim) <- c('opponent','idx','similarity','oppo.farm','itself.farm')
head(msim)

simi <- msim
simi$oppo.farm <- as.factor(simi$oppo.farm)
simi$itself.farm <- as.factor(simi$itself.farm)
levels(simi$oppo.farm) <- c('MN1','MN2','MN6','MN7','MN8','NY','PA','VT')
levels(simi$itself.farm) <- c('MN1','MN2','MN6','MN7','MN8','NY','PA','VT')
simi$flabels <- paste0(simi$itself.farm,'_',simi$oppo.farm)
simi.list <- SplitData(simi,'itself.farm')

# check how many samples in each farm
# nrow(simi.list[[1]])/525 = 38 number of samples in that farm 

qcutoff = 0.99 # 0.98 only works for MN farms
high.sim.list <- list()

for (f in 1:length(sim.list)) {

	fsim <- simi.list[[f]]
	fsim$flabels <- as.factor(fsim$flabels)	
	# density plot to look at the overall distribution
        dens.title <- paste0('Density plot to examine distribution of similarities from pair-wise comparisons ',fsim$itself.farm[1])
	densp <- ggplot(fsim,aes(x=similarity, fill=flabels)) + geom_density(alpha=0.45) + facet_wrap(~ flabels,ncol = 1) +
		 labs(title = dens.title, x = 'Comparisons')+ guides(fill=FALSE)
	densp.title <- paste0(fpath,gsub(" ", "_",dens.title, fixed = TRUE),'.pdf')
        ggsave(densp,file = densp.title,scale = 0.8) # scale <1 larger, >1 smaller
        print(densp)

	q <- quantile(fsim$similarity, probs = c(qcutoff, 1))
	fq <- filter(fsim,(similarity > q[1] & similarity < 0))
	# dot plot to look at the highly similar samples
	dot.title <- paste0('Dotplot to examine highly similar pairs at top ',paste0(100*(1-qcutoff),' percent '),fq$itself.farm[1])
	dotp <- ggplot(fq, aes(x=flabels, y=similarity, fill = flabels)) + geom_dotplot(binaxis='y', stackdir='center', binwidth = 1.2, dotsize=0.1) +
		labs(title = dot.title, x = 'Comparisons')+ guides(fill=FALSE)
	dotp.title <- paste0(fpath,gsub(" ", "_",dot.title, fixed = TRUE),'.pdf')
        ggsave(dotp,file = dotp.title,scale = 0.8) # scale <1 larger, >1 smaller
        print(dotp)
	
	high.sim.list[[f]] <- fq
}


# histogram like:
        dotp <- ggplot(fq, aes(x=flabels, y=similarity, fill = flabels)) + geom_dotplot(binaxis='y', method="histodot", binwidth = 1.2, dotsize=0.3) +
                labs(title = dot.title, x = 'Comparisons')+ guides(fill=FALSE)

















	
