
library(plyr)

PlotHistByVar <- function(total,var,ncols,xcor,ycor,plot.title) { 
                tot <- as.data.frame(cbind(var,total))
                colnames(tot) <- c('Types','Total_SNPs')
                tot$Type <- factor(tot$Types,labels = levels(var))

                summ <- ddply(tot, 'Types', summarise, N=paste('N=', length(Total_SNPs)), Mean= paste('Mean=', round(mean(Total_SNPs),0)))
                summ$Type <- factor(summ$Types,labels = levels(var))

                fw <- ggplot(tot, aes(x = Total_SNPs, fill = Type, color = Type)) +
                        geom_histogram(bins = 50) +
                        ggtitle(plot.title) +
                        facet_wrap(~ Type, ncol = ncols) +
                        geom_text(data = summ,aes(x =  xcor, y = ycor, label = paste(N,Mean, sep = ', ')), inherit.aes = FALSE)
                trim.titl <- paste0(fpath,gsub(" ", "_",plot.title, fixed = TRUE),'.pdf')
                ggsave(fw,file = trim.titl,scale = 1) # scale <1 larger, >1 smaller
                print(fw)
}

# how to position the text at upper right corner for each facet?
# geom_text(data = nobs,aes(x =  -Inf, y = Inf, label = n), hjust = 1, vjust = 1, inherit.aes = FALSE)

