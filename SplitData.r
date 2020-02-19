
SplitData <- function(gen,var.name){
        var <- gen[,which(colnames(gen) == var.name)]
        fact <- levels(as.factor(var))
	split.dat <- list()
        for (d in 1:length(fact)){
                split.dat[[d]] <- filter(gen,var == fact[d])
        }
        return(split.dat)
}

