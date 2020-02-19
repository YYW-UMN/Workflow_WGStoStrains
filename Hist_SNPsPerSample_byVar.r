Hist_SNPsPerSample_byVar <- function(gen,var.name,snp.type = c('total','het','hom'),pname){
	var <- gen[,which(colnames(gen) == var.name)]
        totals <- SumSNPs(gen) # a list of 3 lists: total, het and hom
	if (snp.type == 'total') {
#		tota <- as.data.frame(cbind(var,totals[[1]]))
		titl <- paste0('Histogram Of Total SNPs Per Sample',pname)
		PlotHistByVar(totals[[1]],var, ncols = 3, xcor = 50, ycor = 25, plot.title = titl)
	} else { 
		if (snp.type == 'het') {
	                titl <- paste0('Histogram Of Heterozygous SNPs Per Sample',pname)
			PlotHistByVar(totals[[2]], var, ncols = 3, xcor = 40, ycor = 40, plot.title = titl)
		} else { # hom
			titl <- paste0('Histogram Of Homozygous SNPs Per Sample',pname)
			PlotHistByVar(totals[[3]], var, ncols = 3, xcor = 25, ycor = 45, plot.title = titl)
		}
	}
}





