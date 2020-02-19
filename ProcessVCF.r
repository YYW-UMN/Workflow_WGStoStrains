library(stringr)
library(seqinr)
library(Rsamtools)
library(plyr)
library(dplyr)

ProcessVCF <- function(vfilepath) {
# 1. Import file 
  InputVCF <- vfilepath
#  InputVCF <- '~/mixed/original/Map_K10_VBC_CLI.whitelist.noindel.vcf'
  input <- read.table(file = InputVCF)


  cat('The VCF file imported contains',nrow(input),'SNPs from',ncol(input),'animals.')
  # Output: [1] 796 534 
  ## Note: 796 rows = 796 SNP positions
  ##       534 columns = 9 header columns + 525 sample columns
  ###### input[1,1:10]
  ## Purpose: First nine columns are metadata, data starts at column No.10 
  # Output:  V1 (CHROM)       V2(POS)  V3(ID) V4(REF) V5(ALT) V6(QUAL)     V7(Filter)  V8(INFO)   V9(Genotype Field)       V10(!!!Samples starts at this column)
  #       1  Map_K10_VBC_CLI  356      .      C       T       2.23073e+03  .           see below  GT:DP:AD:RO:QR:AO:QA:GL  0/0:71:71,0:71:2772:0:0:0,-21.3731,-249.62
  # INFO: AB;ABP;AC;AF;AN;AO;CIGAR;DP;DPB;DPRA;EPP;EPPR;GTI;LEN;MEANALT;MQM;MQMR;NS;NUMALT;ODDS;PAIRED;PAIREDR;PAO;PQA;PQR;PRO;QA;QR;RO;RPL;RPP;RPPR;RPR;RUN;SAF;SAP;SAR;SRF;SRP;SRR;TYPE
  pos <- input[,2]
  ref <- input[,4]
  alt <- input[,5]

  ###### str_count(input[1,9], ":")
  # Output: 7 
  ## Purpose: what's stored in input column No.9 (Genotype Field)? 
  ## Note: 7 ":" seperates 8 fields 
  ##       GT(Genotype): 0/0,0/1,
  ##       DP(Read Depth):  
  ##       AD(Number of Observation for each allele):
  ##       RO(Reference allele observation count): 
  ##       QR(Sum of quality of the reference observations): 
  ##       AO(Alternate allele observation count):
  ##       QA(Sum of quality the alternate observations): 
  ##       GL(Genotype likelihood).
  ### Function: str_count() count the occurence of ":" in column 9 of row 1

# 2. Create a vcf file in R from input file
  vcf <- as.matrix(read.table(file = InputVCF,comment.char="",sep="\n"))
  ###### dim(vcf)
  # Output: [1] 858   1 
  ## Note: 858 rows = 62 header rows + 796 SNP rows 
  vcf_header_end_row <- which(grepl("#CHROM",vcf)==TRUE)
  # Output: [62,] 
  ## Note: the header always ends with "#CHROM"

  names_header_line <- unlist(strsplit(vcf[vcf_header_end_row],"\t"))
  # Output: "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMAP129\tMAP381
  ## Note: vcf[vcf_header_end_row] looks like this: "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMAP129\tMAP381
  ##       it's a line with 9 fixed mandatory header line names, and followed by sample names.
  ##       Each name is seperated by '\t' 
  ##       The 9 header names are "#CHROM""POS""ID""REF""ALT""QUAL""FILTER""INFO""FORMAT"
  ### Function: strsplit(string, split by ""),unlist: flatten list [[1]]-> [1]

  ###### length(names_header_line)
  # Output: 534
  ## Note: 534 = 9 headers + 525 samples 

  format <- which(names_header_line=='FORMAT')
  # Output: 9
  ## Purpose: where is the end of header names
  ## Note: In Freebayes, header names end with "FORMAT" 
  
  sample_names <- names_header_line[-c(1:format)]
  # Output: MAP129 - 775

# 3.a Summarize SNPs metrics: the 8th column of input file. Sample data is excluded
  tmp.snp <- data.frame(do.call('rbind', strsplit(as.character(input[,8]),';',fixed=TRUE)))
  snp <- data.frame(input[,2])
  for (s in 1:ncol(tmp.snp)) {
    t <- data.frame(do.call('rbind', strsplit(as.character(tmp.snp[,s]),'=',fixed=TRUE)))
    snp <- cbind(snp,t[,2])
    colnames(snp)[s+1] <- as.character(t[1,1])
  }
  colnames(snp)[1] <- 'position'
  #return this matrix for further analysis 

# 3.b Summarize sample data using eight data.frames ("GT" "DP" "AD" "RO" "QR" "AO" "QA" "GL")
  info <- unlist(str_split(input[1,format], ":"))
  # Output: [1] "GT" "DP" "AD" "RO" "QR" "AO" "QA" "GL"
  ## Purpose:  extract info headers from the Genotype field V9

  n <- length(info) 
  # Output: 8 
  ## Purpose: know how many entries in the GT (genotype) field

# 3.b Step 1: seperate sample data from input
  samples <- input[,c((format+1):ncol(input))]
  ####### dim(samples)
  ## Notes: rows = 796(snps), columns = 525(samples)
  ## Purpose:  extract data of all samples from V10 to the last column

# 3.b Step 2: initiate data.frames for results
  # (1) GT Genotype
  gt <- matrix(0,nrow(samples),ncol(samples))
  # (2) DP Read depth at this position for this sample 
  dp <- matrix(0,nrow(samples),ncol(samples))
  # (3) AD Number of observation for each allele 
  ad <- matrix(0,nrow(samples),ncol(samples))
  # (4) RO Reference allele observation count
  ro <- matrix(0,nrow(samples),ncol(samples))
  # (6) AO Alternate allele observation count 
  ao <- matrix(0,nrow(samples),ncol(samples)) 
  # (5) QR Sum of quality of the reference observation
#  qr <- matrix(0,nrow(samples),ncol(samples))
  # (7) QA Sum of quality of the alternate observation
#  qa <- matrix(0,nrow(samples),ncol(samples)) 
  # (8) GL Genotype Likelihood, 
  # log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy 
#  gl <- matrix(0,nrow(samples),ncol(samples))

# 3.c Step 3: Fill in data
  ### Function: grepl(pattern, string to search) when fixed = TRUE, literally pattern; otherwise regex pattern (for example, '.' is a wild card)
  for (i in 1:nrow(samples)){
    for (j in 1:ncol(samples)){
      sample <- unlist(str_split(samples[i,j], ":"))
      gt[i,j] <- sample[1]
      if (grepl('.',sample[2], fixed = TRUE) == TRUE) {
        dp[i,j] <- 999
        } else { dp[i,j] <- as.numeric(sample[2]) }
      if (grepl('.',sample[4], fixed = TRUE) == TRUE) {
        ro[i,j] <- 999
        } else { ro[i,j] <- as.numeric(sample[4]) }
#      if (grepl('.',sample[5], fixed = TRUE) == TRUE) {
#        qr[i,j] <- 999
#        } else { qr[i,j] <- as.numeric(sample[5]) }
      ad[i,j] <- sample[3] # multiple values seperated by ','
      ao[i,j] <- sample[6] # multiple count for each alt allele, values  seperated by ','
#      qa[i,j] <- sample[7] # multiple quality for each alt allele, values  seperated by ','
#      gl[i,j] <- sample[8] # multiple likelihood for each alt allele, values  seperated by ','
    } 
  }
  cat('writing 8 files for GT, DP, AD, RO, AO, QR, QA, GL...')
  ###### dim(gt)
  # Output: [1] rows = 796 SNPs, column =  525 samples

# 3.c Step 4: Need to do: In ad, ao, qa, gl, each row is a string seperated by ','. 
  fpath <- '~/mixed/output/'
  write.csv(file = paste0(fpath,'gt',nrow(gt),'.csv'),gt) 
  write.csv(file = paste0(fpath,'DP',nrow(dp),'.csv'),dp)
  write.csv(file = paste0(fpath,'AD',nrow(ad),'.csv'),ad)
  write.csv(file = paste0(fpath,'RO',nrow(ro),'.csv'),ro)
  write.csv(file = paste0(fpath,'AO',nrow(ao),'.csv'),ao) 
#  write.csv(file = paste(fpath,'QR.CSV',sep = ''),qr) 
#  write.csv(file = paste(fpath,'QA.csv',sep = ''),qa)
#  write.csv(file = paste(fpath,'GL.csv',sep = ''),gl)
  write.csv(file = paste(fpath,'sampleIDs',nrow(gt),'.csv',sep = ''),sample_names)
  write.csv(file = paste(fpath,'positions',nrow(gt),'.csv',sep = ''),pos)
  cat('Output files generated: Rows are SNPs. Columns are Samples.')
#  return(samples,pos)
	

