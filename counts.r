

# sam X998 X999 X1 X2 X3 X4 X10 X12 X20 X23 good bad
csam <- read.csv(paste(fpath,'GT_filtered_Samples_count.csv',sep = ''),header = TRUE)

# snp X998 X999  X1 X2 X3 X4 X10 X12 X20 X23 good bad hom het
csnp <- read.csv(paste(fpath,'GT_filtered_SNP_count.csv', sep = ''),header = TRUE)

# at least three animals
rows <- which(csnp$good >= 3)

do the summary table to double check the results

meta <- read.csv(paste(mpath,'MetaDATAJohnes2.csv',sep = ''),header = TRUE); meta <- meta[,-1]

ny <- filter(meta, State == "NY")
pa <- filter(meta, State == "PA")
vt <- filter(meta, State == "VT")
mn <- filter(meta, State == "MN")

state <- cbind(meta$State,csam)
colnames(state)[1] <- c('State')

# not filtered 
gt <- read.csv(paste(fpath,'Filtered_GT.csv',sep = ''),header = TRUE)
t <- t(gt)
