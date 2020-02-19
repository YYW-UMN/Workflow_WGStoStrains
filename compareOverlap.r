gt <- read.csv('~/mixed/output/GT_2755_sample_by_snp_with_meta_cols.csv', header = TRUE); gt <- gt[,-c(1,2)]
meta.cols <- GetMetaColsIndex(gt)
meta <- gt[,meta.cols]
meta$idx <- seq(1,nrow(meta),1)
g <- gt[,-meta.cols]
ann <- which(meta$Contamination == 1)

overlap <- matrix(0,nrow = nrow(g),ncol = nrow(g))
# loop over entire matrix
#for (i in 1:nrow(overlap)){
#	baseline <- g[i,]
#	for (j in 1:ncol(overlap)){
#		overlap[i,j] <- length(which(baseline == g[j,]))
#        }
#}

# loop over the diagnal and upper triangle
for (i in 1:nrow(overlap)){
        baseline <- g[i,]
        for (j in i:ncol(overlap)){
		snps <- which(baseline>0)
                overlap[i,j] <- length(which(baseline[,snps] == g[j,snps]))
        }
}

op <- overlap
write.csv(op,'overlap.csv')
# copy the upper trianlge to the lower triangle 
op[lower.tri(op)] <- t(op)[lower.tri(op)]
isSymmetric(op)

mop <- melt(op)
colnames(mop) <- c('idx','itself','snp.overlap')
mop <- join(mop,meta[,c('idx','Farm')],by = 'idx')

colnames(mop) <- c('opponent','idx','snp.overlap','oppo.farm')
mop <- join(mop,meta[,c('idx','Farm')],by = 'idx')

# finalize the column names
colnames(mop) <- c('opponent','idx','snp.overlap','oppo.farm','itself.farm')
head(mop)

movlp <- mop
movlp$oppo.farm <- as.factor(movlp$oppo.farm)
movlp$itself.farm <- as.factor(movlp$itself.farm)
levels(movlp$oppo.farm) <- c('MN1','MN2','MN6','MN7','MN8','NY','PA','VT')
levels(movlp$itself.farm) <- c('MN1','MN2','MN6','MN7','MN8','NY','PA','VT')
movlp$flabels <- paste0(movlp$itself.farm,'_',movlp$oppo.farm)
movlp.list <- SplitData(movlp,'itself.farm')

colnames(movlp)[3] <- c('similarity')

totals <- filter(movlp,opponent==idx)
tot <- rep(0,nrow(movlp))
for (i in 1:nrow(movlp)){
	lookfor <- movlp[i,2]
	tot[i] <- totals$similarity[which(totals$idx == lookfor)]
}

movlp <- cbind(movlp,tot)
write.csv(movlp,'movlp.csv')
movlp$similarity2 <- movlp$similarity/movlp$tot
movlp <- filter(movlp, opponent!=idx)

summary(movlp$similarity2)

qcutoff = 0.99 # 0.98 only works for MN farms
hi.o <- filter(movlp,similarity2 > qcutoff)
dim(hi.o)

movlp.list <- SplitData(hi.o,'itself.farm')

hi.ann <- join(hi.o,meta[,c('idx','Contamination')],by = 'idx')
head(hi.ann[order(-hi.ann$Contamination), ])

nodes <- meta[,c(13,4)]

high.movlp.list <- list()
for (f in 1:length(movlp.list)) {
	fsim <- movlp.list[[f]]
	fsim$flabels <- as.factor(fsim$flabels)
#        q <- quantile(fsim$similarity, probs = c(qcutoff, 1))
        fq <- filter(fsim,similarity2 > qcutoff)
	
        # dot plot to look at the highly similar samples
        hist.title <- paste0('Histogram of samples sharing at least ',paste0(qcutoff*100,' percent of SNPs with each other in '),fq$itself.farm[f])
	
ggplot(fq, aes(x = flabels, fill = flabels)) + geom_bar(stat = 'count',width = 0.7)+ labs(title = hist.title, x = 'State labels from highly similar pairs')+ guides(fill=FALSE)

dd <-cbind(fq$idx, fq$opponent)
gg <- graph.data.frame(dd, directed=FALSE)

nodes <- union(unique(fq$idx),unique(fq$opponent))
nmt <- filter(meta,idx %in% nodes)
V(gg)$nmt[,4]
plot(gg, vertex.size = 20)

histp <- ggplot(fq, aes(x=similarity2, fill = flabels)) + geom_histogram(binwidth=2) + facet_grid(~ flabels)+ labs(title = hist.title, x = 'Proportion of identical SNPs')+ guides(fill=FALSE)	
        histp.title <- paste0(fpath,gsub(" ", "_",hist.title, fixed = TRUE),'.pdf')
        ggsave(histp,file = histp.title,scale = 0.8) # scale <1 larger, >1 smaller
        print(histp)
        high.movlp.list[[f]] <- fq
}


        dotp <- ggplot(fq, aes(x=flabels, y=similarity2*100, fill = flabels)) + geom_dotplot(binaxis='y', stackdir='center', binwidth = 1.2, dotsize=0.05) +labs(title = dot.title, x = 'Comparisons')+ guides(fill=FALSE)




#### Look only homozygous SNPs 

gt$idx <- seq(1,nrow(gt),1)
cornel <- filter(gt, State != 'MN')

ct <- CountSNPs(cornel) # count 1,2,3,4,10,20,12,23,998,999,good,bad,hom,het,mix for each SNP
length(which(ct['mixed',] == 0))
snpss <- which((ct['het',] == 0)&(ct['bad',] == 0))
cornell <- cornel[,snpss+1] # because SAMPLE ID column is not in ct
cornell <- cornell[,colSums(cornell)>0]
dim(cornell) # 419 91


overlap <- matrix(0,nrow = nrow(cornell),ncol = nrow(cornell))
for (i in 1:nrow(overlap)){
        baseline <- cornell[i,]
        for (j in i:ncol(overlap)){
                overlap[i,j] <- length(which(baseline == cornell[j,]))
        }
}

op <- overlap
op[lower.tri(op)] <- t(op)[lower.tri(op)]
isSymmetric(op)

library(reshape)
mop <- melt(op)
colnames(mop) <- c('idx','itself','snp.overlap')
mop <- join(mop,meta[,c('idx','Farm')],by = 'idx')

colnames(mop) <- c('opponent','idx','snp.overlap','oppo.farm')
mop <- join(mop,meta[,c('idx','Farm')],by = 'idx')

# finalize the column names
colnames(mop) <- c('opponent','idx','snp.overlap','oppo.farm','itself.farm')
head(mop)
 

movlp <- mop
movlp$oppo.farm <- as.factor(movlp$oppo.farm)
movlp$itself.farm <- as.factor(movlp$itself.farm)
levels(movlp$oppo.farm) <- c('NY','PA','VT')
levels(movlp$itself.farm) <- c('NY','PA','VT')
movlp$flabels <- paste0(movlp$itself.farm,'_',movlp$oppo.farm)
#movlp.list <- SplitData(movlp,'itself.farm')

goto <- which(movlp$idx == movlp$opponent)
identical.n <- movlp$snp.overlap[goto]
# attentiaon: the maximum number of overlap is 91 which is the total number of snps, but this number may change for other datasets
overlap.perct <- round(movlp$snp.overlap/identical.n[1],3)
movlp <- cbind(movlp,overlap.perct)
head(movlp)

#value cutoff
#vcutoff = round(ncol(cornell)*qcutoff) # 0.98 only works for MN farms
#hi.o <- filter(movlp,snp.overlap > vcutoff)

# percentage cutoff
qcutoff = 0.980 # 0.98 only works for MN farms
hi.o <- filter(movlp,overlap.perct > qcutoff)

#movlp.list <- SplitData(hi.o,'itself.farm')
hi.o <- filter(hi.o, idx != opponent)

meta.cornel <- filter(meta, State != 'MN')
hi.ann <- join(hi.o,meta.cornel[,c('idx','Contamination')],by = 'idx')


head(hi.ann[order(-hi.ann$Contamination), ])
sorted.ann <- hi.ann[order(-hi.ann$Contamination), ]

unique(hi.o$idx)
contam <- filter(meta,Contamination == 1) 

sorted.ann.contam <- filter(sorted.ann, Contamination == 1)
sorted.ann.contam

sorted.ann.list <- SplitData(sorted.ann,'itself.farm')

high.movlp.list <- list()
for (f in 1:length(movlp.list)) {
        fann <- sorted.ann.list[[f]]
        fann$flabels <- as.factor(fann$flabels)
#        q <- quantile(fsim$similarity, probs = c(qcutoff, 1))
#        fq <- filter(fsim,similarity > vcutoff)
        # dot plot to look at the highly similar samples
        dot.title <- paste0('Dot plot of samples sharing at least ',qcutoff*100,' percent of SNPs with each other in ',fann$itself.farm[f])
        histp <- ggplot(fann, aes(x=overlap.perct*100, fill = flabels)) + geom_histogram(binwidth=2) + facet_grid(~ flabels)+ labs(title = hist.title, x = 'Identical SNPs')+ guides(fill=FALSE)

        dotp <- ggplot(fann, aes(x=flabels, y=overlap.perct*100, fill = flabels)) + geom_dotplot(binaxis='y', stackdir='center', binwidth = 1.2, dotsize=0.05) +
                labs(title = dot.title, x = 'Comparisons')+ guides(fill=FALSE)
	




        histp.title <- paste0(fpath,gsub(" ", "_",hist.title, fixed = TRUE),'.pdf')
        ggsave(histp,file = histp.title,scale = 0.8) # scale <1 larger, >1 smaller
        print(histp)
        high.movlp.list[[f]] <- fq
}


hi.ann <- join(hi.o,meta[,c('idx','Contamination')],by = 'idx')







