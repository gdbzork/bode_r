loadCDT <- function(fn) {
  print(fn)
  m <- read.table(fn,header=T)
  return(m)
}

cdt2density <- function(cdt,tag) {
  m <- cdt[,4:dim(cdt)[2]]
  d <- apply(m,c(2),mean)
  l <- length(d)
  l2 <- l / 2
  df <- data.frame((-l2+1):l2,d,tag)
  colnames(df) = c('col','mean','tag')
  return(df)
}

library(ggplot2)
fn = commandArgs(T)[1]
out = commandArgs(T)[2]

files = read.table(fn,header=F,stringsAsFactors=F)
colnames(files) = c('tag','fn')
df <- data.frame()
for (i in 1:dim(files)[1]) {
  m <- loadCDT(files[i,2])
  d <- cdt2density(m,files[i,1])
  df <- rbind(df,d)
}
pdf(out)
qplot(col,mean,data=df,color=tag,geom="line")
dev.off()
