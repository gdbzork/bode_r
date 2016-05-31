biasForFile <- function(fn) {
  df <- read.table(fn,header=F)
  colnames(df) <- c('chrom','left','right','name','score','bias','top','bot')
  reads <- sum(df$top) + sum(df$bot)
  bdf <- df[df$bias >= 10,]
  biased <- sum(bdf$top) + sum(bdf$bot)
  return(biased/reads)
}

fns <- commandArgs(T)
outFN <- fns[1]
fns <- fns[2:length(fns)]
results <- sapply(fns,biasForFile)
outDF <- data.frame(results,fns,stringsAsFactors=F)
write.table(outDF,file=outFN,row.names=F,col.names=F,quote=F)
