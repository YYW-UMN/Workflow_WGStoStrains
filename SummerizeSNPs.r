
farms.geno

SummerizeSNPs <- function(farms.geno) {

	for (f in (1:length(farms.geno))){
		geno <- farms.geno[[f]] 
                counts <- matrix(0,ncol = 4,nrow = ncol(geno))
                for (s in (1:ncol(geno))){ 
                        counts[s,1] <- sum(geno[,s] == 1)+ sum(geno[,s] == 2)+ sum(geno[,s] == 3)+ sum(geno[,s] == 4) # hom 
                        counts[s,2] <- sum(geno[,s] == 10)+ sum(geno[,s] == 20)+ sum(geno[,s] == 12)+ sum(geno[,s] == 23) # het
                        counts[s,3] <- sum(geno[,s] == 999) + sum(geno[,s] == 998) # missing value or low coverage
                        counts[s,4] <- sum(geno[,s] != 0) # all SNPs
                }

	colnames(counts) <- c('hom','het','low','all') # sum(counts[,4] < round(animals*percent))
	hom <- list(); het <- list(); mixed <- list(); rare <- list()
	snps <- rep(0,nrow(counts))
	for (s in (1:nrow(counts))) {
		if (counts[s,1] != 0 & counts[s,2] == 0){
			snps[s] <- 0 # hom
		} else {
			if (counts[s,1] == 0 & counts[s,2] != 0){
				snps[s] <- 1 # het
			} else {
				if (counts[s,1] != 0 & counts[s,2] != 0) {
					snps[s] <- 10 # mixed
				} else {
					snps[s] <- 999 # low
				}
			}
		}
	}

	farm <- farms.geno
	pos <- c()
	for (n in (1:length(colnames(farms.geno)))){
		p <- as.numeric(unlist(strsplit(colnames(farms.geno)[n],"X"))[2])
		pos <- append(pos,p)
	}
	farm$Date <- mapvalues(farm$Date,'unknown',999)
	snp.count <- matrix(0, nrow = length(pos),ncol = length(unique(farm$Date)))
	#sample.count <- matrix(0, nrow = length(pos),ncol = length(unique(farm$Date)))

	colnames(snp.count) <- sort(unique(farm$Date))
	rownames(snp.count) <- colnames(geno)

        sample.count <- matrix(0, nrow = length(pos),ncol = length(unique(farm$Date)))
        colnames(sample.count) <- sort(unique(farm$Date))
        rownames(sample.count) <- colnames(geno)
	
	snp.list <- grep('X',colnames(farm))
	date <- grep('Date',colnames(farm))
	#snp.count <- cbind(pos,snps,snp.count)
	#years <- as.numeric(colnames(snp.count))
	dat <- farm[,c(date,start:ncol(farm))]
	

	for (snp in 1:nrow(sample.count)){ # rows are snps
		
			sample.count[snp,] <- aggregate(dat[,cols],by = list(Category = dat[,'Date']), FUN = length)[,2]
		#non.zero.rows <- which(dat[,snp+1] > 0 & dat[,snp+1] < 998) 		
	}










	









