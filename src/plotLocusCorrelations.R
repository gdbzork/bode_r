loadCounts <- function(fn) {
  raw = read.table(fn,header=FALSE,as.is=TRUE,sep="\t")
  merge = raw$V1 + raw$V2 + raw$V3
  return(merge)
}

fn2tag <- function(fns) {
  flds = strsplit(fns,"_")
  code = sapply(flds,function(x) x[1])
  tissue = sapply(flds,function(x) x[3])
  tag = sprintf("%s_%s",code,tissue)
  return(tag)
}

fn2meta <- function(fns) {
  flds = strsplit(fns,"_")
  code = sapply(flds,function(x) x[1])
  tissue = sapply(flds,function(x) x[3])
  species = sapply(flds,function(x) substr(x[5],1,3) )
  mat = cbind(code,tissue,species)
  colnames(mat) = c('code','tissue','species')
  return(data.frame(mat,stringsAsFactors=FALSE))
}

lists2matrix <- function(lset,cnames) {
  uni = vector()
  for (s in lset) {
    uni = union(uni,s)
  }
  mat = cbind()
  for (s in lset) {
    mat = cbind(mat,uni %in% s)
  }
  colnames(mat) = cnames
  rownames(mat) = uni
  return(mat)
}

fn2taglist <- function(fns) {
  flds = strsplit(fns,"[_.]")
  tissue = sapply(flds,function(x) x[3])
  return(tissue)
}
fn2code <- function(fns) {
  flds = strsplit(fns,"[_.]")
  code = sapply(flds,function(x) x[1])
  return(code)
}

relevantReads <- function(fns) {
  stuff = rbind()
  keys = fn2code(fns)
  for (i in 1:length(fns)) {
    r = loadCounts(fns[i])
    co = keys[i]
    s = sum(r)
    stuff = rbind(stuff,c(co,s))
  }
  stuff2 = data.frame(stuff[,1],as.numeric(stuff[,2]),stringsAsFactors=FALSE)
  colnames(stuff2) = c('code','trna')
  return(stuff2)
}
  
loadSamples <- function(fns) {
  tags = fn2taglist(fns)
  res = list()
  for (i in 1:length(fns)) {
    x = read.table(fns[i],header=FALSE,as.is=TRUE)
    res[tags[i]] = x
  }
  return(res)
}

fns = commandArgs(TRUE)

tags = fn2tag(fns)
outfn = sprintf("%s_%s.pdf",tags[1],tags[2])
raw = cbind()
for (fn in fns) {
  raw = cbind(raw,loadCounts(fn))
}
keep = apply(raw,c(1),sum) > 2
print(sum(keep))
filt = raw[keep,]

pdf(outfn)
plot(filt[,1],filt[,2],xlab=tags[1],ylab=tags[2],main="Mouse Tissue Plot")
dev.off()
