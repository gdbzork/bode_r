# Make a table of loci which are not considered transcribed, showing how many loci have zero, 1, 2, ...
# reads covering them (and their flanks).
#
# Inputs:
#   1) unfiltered read counts, for each library
#   2) pol3_library_metadata.Rdata
#   3) lists of "clean" loci and "transcribed" loci
#
# Outputs:
#   1) writes table to "results.txt"

code2fn <- function(code,filelist) {
  splitlist = strsplit(filelist,"_")
  found = sapply(splitlist,function(x) x[1] == code)
  return(filelist[found])
}

loadCounts <- function(fn,libdata) {
  tbl = read.table(fn,header=F,stringsAsFactors=F)
  colnames(tbl) = c('upstream','middle','downstream','chrom','left','right','name','score','strand')
  code = strsplit(fn,"_")[[1]][1]
  species = libdata[libdata$code==code,]$species
  cleanfn = paste("trna",species,"clean.txt",sep="_")
  realfn = paste("trna",species,"names.txt",sep="_")
  cleanList = read.table(cleanfn,header=F,stringsAsFactors=F)
  realList = read.table(realfn,header=F,stringsAsFactors=F)
  clean = tbl[tbl$name %in% cleanList$V1,]
  unreal = clean[!(clean$name %in% realList$V1),]
  unreal$total = unreal$upstream + unreal$middle + unreal$downstream;
  return(unreal)
}

sfunc <- function(a,b,n) {
  a = a+1
  if (a > n) {
    a <- n
  }
  b[a] <- b[a] + 1
  return(a)
}

summarize <- function(counts,n) {
  bins = integer(n)
#  sapply(counts$total,sfunc,bins,n)
  for (i in counts$total) {
    j = i+1
    if (j > n) {
      j <- n
    }
    bins[j] <- bins[j] + 1
  }
  return(bins)
}

filelist <- list.files()
load('pol3_library_metadata.Rdata')
relevant <- pol3_library_meta[pol3_library_meta$inPaper,]
bins <- 11
labels <- paste((1:bins) - 1,sep="")
outtbl <- cbind(t(labels),"code","filename")
for (code in relevant$code) {
  fn <- code2fn(code,filelist)
  counts <- loadCounts(fn,relevant)
  s <- summarize(counts,bins)
  r <- cbind(t(paste(s,sep="")),code,fn)
  outtbl <- rbind(outtbl,r)
}
write.table(outtbl,file="results.txt",col.names=F,row.names=F,quote=F)
