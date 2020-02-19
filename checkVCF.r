
library(VariantAnnotation)
library(BSgenome.MapK10.SRA.SRR060191.3)

vcf_files <- c(
	'~/mixed/original/Map_K10_VBC_CLI.0.2000000.vcf.gz',
	'~/mixed/original/Map_K10_VBC_CLI.2000000.4000000.vcf.gz',
	'~/mixed/original/Map_K10_VBC_CLI.4000000.4832589.vcf.gz'
)

vcf_list <- lapply(vcf_files, function(vcf_file){
	cat(sprintf('reading %s\n', vcf_file))
	vcf <- readVcf(vcf_file, 'MapK10')
	#is_frequent_snp <- rowSums(geno(vcf)$GT != '0/0') >= 1
	is_frequent_snp2 <- rowSums(geno(vcf)$GT != '0/0' & !is.na(geno(vcf)$DP) & geno(vcf)$DP > 20) >= 1
	is_on_forward_strand <- sapply(info(vcf)$SAF, function(x) all(x > 0))
	is_on_reverse_strand <- sapply(info(vcf)$SAR, function(x) all(x > 0))
	good <- is_frequent_snp2 & is_on_forward_strand & is_on_reverse_strand
	vcf[good]
})


v <- Reduce('rbind', vcf_list)
output_file <- '~/mixed/original/Map_K10_VBC_CLI.rds'
saveRDS(v, output_file)

output_file <- '~/mixed/original/Map_K10_VBC_CLI.vcf'
writeVcf(v, output_file)

bed_file <- '~/mixed/original/exclude_all.bed'
bed <- read.table(bed_file, sep = '\t', header = TRUE)
bed <- GRanges(seqnames = bed[, 1], range = IRanges(bed[, 2], bed[, 3]))
v <- v[!v %over% bed]
output_file <- '~/mixed/original/Map_K10_VBC_CLI.whitelist.rds'
saveRDS(v, output_file)
output_file <- '~/mixed/original/Map_K10_VBC_CLI.whitelist.vcf'
writeVcf(v, output_file)

is_indel <- width(rowRanges(v)$REF) > 1 | sapply(rowRanges(v)$ALT, function(alt) sum(width(alt) == 1) == 0)
output_file <- '~/mixed/original/Map_K10_VBC_CLI.whitelist.noindel.rds'
saveRDS(v, output_file)
output_file <- '~/mixed/original/Map_K10_VBC_CLI.whitelist.noindel.vcf'
writeVcf(v, output_file)


is_frequent_snp <- rowSums(geno(v)$GT != '0/0' & !is.na(geno(v)$DP) & geno(v)$DP > 20) >= 2

# top: number of samples; bottom: SNPs  
xx <- rowSums(geno(v)$GT != '0/0')
table(xx)

# top: number of samples with '.', bottom: SNPs 
a <- rowSums(geno(v)$GT == '.')
v <- readRDS('~/mixed/original/Map_K10_VBC_CLI.whitelist.noindel.rds')















