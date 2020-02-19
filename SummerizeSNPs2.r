SummarizeSNPs2 <- function(cdat){
        if (typeof(cdat) == "list")  { # more than one data
                sum.table <- matrix(0,ncol = length(cdat), nrow = 4)    
                for (i in 1:length(cdat)){
                        c <- cdat[[i]]
                        remov <- which(c[12,] == 0)
                        if (length(remov) > 0) {
				c <- c[,-remov]
			}
                        sum.table[1,i] <- length(which(c[13,] != 0)) # hom
                        sum.table[2,i] <- length(which(c[14,] != 0)) # het
                        sum.table[3,i] <- length(which(c[15,] != 0)) # mix
			sum.table[4,i] <- length(which(c[12,] != 0)) # total
                }
                return(sum.table)
        } else {                       
                remov <- which(cdat[12,] == 0) #(if good == 0, then remove snp from analysis)
                if (length(remov) > 0) {
                	cdat <- cdat[,-remov]
                }
                sum.row <- rep(0,3)                                     
                sum.row[1] <- length(which(cdat[13,] != 0))# hom
                sum.row[2] <- length(which(cdat[14,] != 0))# het
                sum.row[3] <- length(which(cdat[15,] != 0))# mixed
                sum.row[4] <- length(which(cdat[12,] != 0))# total
                return(sum.row)
        }
}       

