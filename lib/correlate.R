library('preprocessCore')
count_cols = c('left','overlap','right','chrom','start','end','name','cove','strand')
summ_cols = c('left','overlap','right','name')

loadCounts <- function(fn) {
  c = read.table(fn,header=FALSE)
  colnames(c) = count_cols
  c
}

loadSumms <- function(fn) {
  c = read.table(fn,header=FALSE)
  colnames(c) = summ_cols
  c
}

tsum <- function(a) {
  a$left+a$overlap+a$right
}

fn2key <- function(name) {
  flds = strsplit(name,"_")[[1]]
  paste(c(flds[1],flds[3],substring(flds[5],1,3)),collapse="_")
}

loadAll <- function(fns) {
  names = sapply(fns,fn2key,USE.NAMES=FALSE)
  objs = sapply(fns,loadSumms,USE.NAMES=FALSE)
  mat = apply(objs,c(2),tsum)
  rownames(mat) = objs[,1]$name
  colnames(mat) = names
  mat
}

tcor <- function(a,b,an,bn) {
  ac = a$left+a$overlap+a$right
  bc = b$left+b$overlap+b$right
  sp = cor(ac,bc,method="spearman")
  pe = cor(ac,bc,method="pearson")
  stxt = sprintf("spearman R=%6.4f  R^2=%6.4f",sp,sp*sp)
  ptxt = sprintf("pearson R=%6.4f  R^2=%6.4f",pe,pe*pe)
  plot(ac,bc,xlab=an,ylab=bn,main="Raw")
  amax = max(ac) * 0.50
  bmax = max(bc) * 0.10
  text(c(amax,amax),c(bmax,bmax/2),pos=4,labels=c(stxt,ptxt))
}

tcorlog <- function(a,b,an,bn) {
  ac = log2(a$left+a$overlap+a$right)
  bc = log2(b$left+b$overlap+b$right)
  sp = cor(ac,bc,method="spearman")
  pe = cor(ac,bc,method="pearson")
  stxt = sprintf("spearman R=%6.4f  R^2=%6.4f",sp,sp*sp)
  ptxt = sprintf("pearson R=%6.4f  R^2=%6.4f",pe,pe*pe)
  plot(ac,bc,xlab=an,ylab=bn,xlim=c(0,max(ac)*1.1),ylim=c(0,max(bc)*1.1),main="Log2")
  amax = max(ac) * 0.50
  bmax = max(bc) * 0.20
  print(amax)
  print(bmax)
  text(c(amax,amax),c(bmax,bmax/2),pos=4,labels=c(stxt,ptxt))
}

tcorn <- function(a,b,an,bn) {
  ac = a$left+a$overlap+a$right
  bc = b$left+b$overlap+b$right
  cb = cbind(ac,bc)
  m = as.matrix(cb)
  mnorm = normalize.quantiles(m,copy=TRUE)
  acnorm = mnorm[,1]
  bcnorm = mnorm[,2]
  sp = cor(acnorm,bcnorm,method="spearman")
  pe = cor(acnorm,bcnorm,method="pearson")
  stxt = sprintf("spearman R=%6.4f  R^2=%6.4f",sp,sp*sp)
  ptxt = sprintf("pearson R=%6.4f  R^2=%6.4f",pe,pe*pe)
  plot(acnorm,bcnorm,xlab=an,ylab=bn,main="Quantile Normalized")
  amax = max(acnorm) * 0.50
  bmax = max(bcnorm) * 0.10
  text(c(amax,amax),c(bmax,bmax/2),pos=4,labels=c(stxt,ptxt))
}

tcorpair <- function(a,b,an,bn) {
  ac = a
  bc = b
  sp = cor(ac,bc,method="spearman")
  pe = cor(ac,bc,method="pearson")
  stxt = sprintf("spearman R=%6.4f  R^2=%6.4f",sp,sp*sp)
  ptxt = sprintf("pearson R=%6.4f  R^2=%6.4f",pe,pe*pe)
  plot(ac,bc,xlab=an,ylab=bn)
  amax = max(ac) * 0.50
  bmax = max(bc) * 0.10
  text(c(amax,amax),c(bmax,bmax/2),pos=4,labels=c(stxt,ptxt))
}


