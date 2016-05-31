loadCDT <- function(fn) {
  print(fn)
  m <- read.table(fn,header=T)
  return(m)
}

cdt2density <- function(cdt1,cdt2,tag) {
  m1 <- cdt1[,4:dim(cdt1)[2]]
  m2 <- cdt2[,4:dim(cdt2)[2]]
  m <- rbind(m1,m2)
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
df <- data.frame()
for (i in 1:dim(files)[1]) {
  m1 <- loadCDT(files[i,2])
  m2 <- loadCDT(files[i,3])
  d <- cdt2density(m1,m2,files[i,1])
  df <- rbind(df,d)
}
pdf(out)
qplot(col,mean,data=df,color=tag,geom="line")
dev.off()
