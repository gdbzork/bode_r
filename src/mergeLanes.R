#!/usr/bin/Rscript

# Quick script to merge K lanes of Solexa data, either per-locus counts or
# grouped by isotype or anticodon.  Use quantile normalization to ensure that
# lanes are weighted equally.  Also sum up the flanks and overlapping regions.

loadSumms <- function(fn) {
  read.table(fn,header=FALSE)
}

tsum <- function(a) {
  a$V1+a$V2+a$V3
}

fn2key <- function(name) {
  flds = strsplit(name,"_")[[1]]
  paste(c(flds[1],flds[3],substring(flds[5],1,3)),collapse="_")
}

loadAll <- function(fns) {
  names = sapply(fns,fn2key,USE.NAMES=FALSE)
  sapply(fns,loadSumms,USE.NAMES=FALSE)
}

collapse <- function(grid) {
  mat = apply(grid,c(2),tsum)
  mat_norm = normalize.quantiles(mat)
  apply(mat_norm,c(1),mean)
}

library(preprocessCore)
fns <- commandArgs(TRUE)
grid <- loadAll(fns)
coll <- collapse(grid)
t <- grid[,1]
t$V1 = NULL
t$V2 = NULL
t$V3 = NULL
x = data.frame(coll,t)
write.table(x,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
