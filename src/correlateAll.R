#!/usr/bin/Rscript

# Quick script to correlate multiple lanes against each other.

loadSumms <- function(fn) {
  read.table(fn,header=FALSE)
}

fn2key <- function(name) {
  flds = strsplit(name,"_")[[1]]
  paste(c(flds[1],flds[2],flds[3]),collapse="_")
}

loadAll <- function(fns) {
  sapply(fns,loadSumms,USE.NAMES=FALSE)
}

library(preprocessCore)
fns <- commandArgs(TRUE)
grid <- loadAll(fns)
names <- sapply(fns,fn2key,USE.NAMES=FALSE)

compnames = c()
spear = c()
pears = c()
for (i in 1:(length(fns)-1)) {
  for (j in (i+1):length(fns)) {
    tag = sprintf("%s_vs_%s",names[i],names[j])
    suff = strsplit(fns[i],"_")[[1]][4]
    compnames = c(compnames,tag)
    plotfn = sprintf("%s_%s.pdf",tag,suff)
    sample1 = grid[,i]
    sample2 = grid[,j]
    spear = c(spear,cor(sample1$V1,sample2$V1,method='spearman'))
    pears = c(pears,cor(sample1$V1,sample2$V1,method='pearson'))
    pdf(plotfn)
    plot(sample1$V1,sample2$V1,xlab=names[i],ylab=names[j],main=tag)
    abline(lsfit(sample1$V1,sample2$V1))
    warnings()
    dev.off()
  }
}
spsq = spear * spear
pesq = pears * pears

out = cbind(compnames,spear,spsq,pears,pesq)
colnames(out) = c('comparison','Spearman','Spearman^2','Pearson','Pearson^2')
write.table(out,col.names=TRUE,quote=FALSE,row.names=FALSE,sep="\t")
