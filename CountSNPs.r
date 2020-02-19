CountSNPs <- function(gen){
	start <- grep('X',colnames(gen))[1]
	end <- tail(grep('X',colnames(gen)), n =1)
        count <- matrix(nrow = 15, ncol = length(colnames(gen)[start:end]),0)
        colnames(count) <- colnames(gen)[start:end]
        rownames(count) <- c(999,998,1,2,3,4,10,20,12,23,'bad','good','hom','het','mixed')
        dat <- gen[,start:end]
        for (i in 1:ncol(dat)) { # each i is a snp
                count[1,i] <- sum(dat[,i] == 999)
                count[2,i] <- sum(dat[,i] == 998)
                count[3,i] <- sum(dat[,i] == 1)
                count[4,i] <- sum(dat[,i] == 2)
                count[5,i] <- sum(dat[,i] == 3)
                count[6,i] <- sum(dat[,i] == 4)
                count[7,i] <- sum(dat[,i] == 10)
                count[8,i] <- sum(dat[,i] == 20)
                count[9,i] <- sum(dat[,i] == 12)
                count[10,i] <- sum(dat[,i] == 23) 
	}

        for (i in 1:ncol(count)) {  
	       	count[11,i] <- sum(count[1:2,i]) # bad
                count[12,i] <- sum(count[3:10,i]) # good, this is also the total number of SNPs
                hom <- sum(count[3:6,i]) # hom
                het <- sum(count[7:10,i]) # het
		if (hom > 0) {
			if (het == 0) {
				count[13,i] <- 1
			} else {
				count[15,i] <- 1
			}
		} else {
			count[14,i] <- 1
        	}
	}
        return(count)
}


