
library(biofiles)
library(VariantAnnotation)

gbkfile <- '/Users/yuanyuanwang/Dropbox/mixed/original/Map_K10_VBC_CLI.gbk'
gb <- biofiles::gbRecord(gbkfile, TRUE)
biofiles::summary(gb)
gb_ranges <- biofiles::ranges(gb, include = c("locus_tag", "protein_id", "product"))

v <- readRDS('~/mixed/original/Map_K10_VBC_CLI.whitelist.noindel.rds')
vr <- rowRanges(v)

seqlevels(gb_ranges) <- "Map_K10_VBC_CLI"
seqnames(gb_ranges) <- "Map_K10_VBC_CLI"
genome(vr) <- "Mycobacterium avium subsp. paratuberculosis strain K-10 [VBC/CLI]" 

mygenes <- findOverlaps(vr, gb_ranges)
