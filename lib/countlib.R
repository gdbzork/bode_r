# A library of R functions to process counts lists, turning them into useful
# R objects and corresponding data files, and generating summaries, etc.

library(preprocessCore)

# a set of aliases, for cases where we've changed TF names along the way,
# or used multiple names for the same thing
aliasSet <- data.frame(key=c('RNAP3','SeC(e)'),val=c('PolIII','SeC'))
rownames(aliasSet) <- aliasSet$key

# A couple of functions to map aliases to core terms, taking vectors or
# lists respectively
mapAlias <- function(words) {
  vals <- as.vector(aliasSet[words,'val'])
  ifelse(is.na(aliasSet[words,'val']),words,vals)
}

mapAliasList <- function(words) {
  sapply(words,mapAlias)
}

# Convert a vector of filenames to a data frame of info gleaned from the
# filename.
filename2annot <- function(filenames) {
  fldset <- strsplit(filenames,"_")
  code <- sapply(fldset,function(x) x[1])
  factor <- mapAliasList(sapply(fldset,function(x) x[2]))
  species <- sapply(fldset,function(x) substring(x[5],1,3))
  tissue <- sapply(fldset,function(x) x[3])
  lane <- sapply(fldset,function(x) x[6])
  tbx <- data.frame(code,factor,tissue,species,lane)
  rownames(tbx)  = filenames
  tbx
}

annot2tags <- function(tags) {
  tags$lane <- NULL
  tags <- unique(tags)
  rownames(tags) <- tags$code
  return(tags)
}

# convert tRNA locus names to isotype and anticodon
locusName2tags <- function(loci) {
  fldset <- strsplit(loci,"_")
  isotype <- mapAliasList(sapply(fldset,function(x) x[1]))
  anticodon <- sprintf("%s_%s",isotype,sapply(fldset,function(x) x[2]))
  df <- data.frame(isotype,anticodon)
  colnames(df) <- c('isotype','anticodon')
  rownames(df) <- loci
  df
}

# load raw counts into a minimal data frame
loadRawCounts <- function(filename) {
  tbl <- read.table(filename,header=FALSE,as.is=c(7))
  df <- data.frame(tbl$V1+tbl$V2+tbl$V3,tbl$V7,stringsAsFactors=FALSE)
  colnames(df) <- c('count','locus')
  df
}

# load lane, return data frame with counts and appropriate annotation
loadLane <- function(filename) {
  tags <- filename2annot(filename)
  raw <- loadRawCounts(filename)
  bio <- locusName2tags(raw$locus)
  rownames(tags) <- NULL
  data.frame(raw,tags,bio)
}

# load data from all filenames in list
#
# Note: must add something to merge multiple lanes of same library
loadAll <- function(fns) {
  tbl <- loadLane(fns[1])
  fncount <- length(fns)
  for (i in 2:fncount) {
    tb <- loadLane(fns[i])
    tbl <- rbind(tbl,tb)
  }
  tbl
}

codes2keys <- function(codes,tags) {
  subtags <- tags[codes,]
  sprintf("%s_%s_%s",subtags$code,subtags$factor,subtags$tissue)
}

raw2counts <- function(tbl,tags,species="mmu",tissue="all") {
  if (tissue == "all") {
    dslice <- tbl[tbl$species==species,c('locus','code','count')]
  } else {
    dslice <- tbl[tbl$species==species & tbl$tissue==tissue,c('locus','code','count')]
  }
  mat <- reshape(dslice,idvar='locus',v.names=c('count'),timevar=c('code'),direction='wide')
  rownames(mat) <- mat$locus
  mat$locus <- NULL
  cn <- colnames(mat)
  codes <- ifelse(substring(cn,1,5)=='count',substring(cn,7),cn)
  colnames(mat) <- codes2keys(codes,tags)
  mat
}

raw2countsRNA <- function(tbl,tags,species="mmu",tissue="all") {
  if (tissue == "all") {
    dslice <- tbl[tbl$species==species,c('locus','code','count')]
  } else {
    dslice <- tbl[tbl$species==species & tbl$tissue==tissue,c('locus','code','count')]
  }
  mat <- reshape(dslice,idvar='locus',v.names=c('count'),timevar=c('code'),direction='wide')
  mat$locus <- NULL
  cn <- colnames(mat)
  codes <- ifelse(substring(cn,1,5)=='count',substring(cn,7),cn)
  colnames(mat) <- codes2keys(codes,tags)
  mat
}

makeNormal <- function(df) {
  dfn <- normalize.quantiles(as.matrix(df))
  rownames(dfn) <- rownames(df)
  colnames(dfn) <- colnames(df)
  dfn
}

aggregateCounts <- function(mat,fld="isotype") {
  tags <- locusName2tags(rownames(mat))
  tbl <- apply(mat,c(1),mean)
  tapply(tbl,INDEX=tags[[fld]],FUN=sum,simplify=TRUE)
}

tabulateTypes <- function(slist,class="isotype") {
  tset <- c()
  for (species in slist) {
    tblname <- sprintf("%s_all",species)
    tags <- locusName2tags(rownames(get(tblname)))
    tset <- c(tset,as.vector(tags[[class]]))
  }
  tset <- unique(sort(tset))
  thing <- c()
  for (species in slist) {
    tblname <- sprintf("%s_all",species)
    tags <- locusName2tags(rownames(get(tblname)))
    tbl <- table(tags[[class]])[tset]
    rownames(tbl) <- tset
    thing <- cbind(thing,tbl)
  }
  colnames(thing) <- slist
  ifelse(is.na(thing),0,thing)
}

mergeSets <- function(tbls,hdrs,relevant,fld="isotype") {
  ds <- c()
  for (t in tbls) {
    hd <- apply(t,c(1),mean)
    tags <- locusName2tags(rownames(t))
    su <- aggregateCounts(as.matrix(hd),fld=fld)
    ds <- cbind(ds,su[relevant])
  }
  colnames(ds) = hdrs
  rownames(ds) = relevant
  ds
}

ackey2isokey <- function(ack) {
  flds <- strsplit(ack,"_")
  sapply(flds,function(x) {x[1]})
}

anticodons2isotypes <- function(acc) {
  keys <- ackey2isokey(rownames(acc))
  lst <- aggregate(acc,list(isotype=keys),FUN=sum)
  rownames(lst) = lst$isotype
  lst$isotype <- NULL
  lst
}

counts2anticodons <- function(mat) {
  keys <- locusName2tags(rownames(mat))
  lst <- aggregate(mat,list(anticodon=keys$anticodon),FUN=sum)
  rownames(lst) <- lst$anticodon
  lst$anticodon <- NULL
  lst
}

makeCorrelationPlot <- function(tbl,years,fn) {
  pdf(fn)
  plot(years,tbl$iso_cor,ylim=c(0,1),xlim=c(0,200),pch=19,xlab="Lineage divergence in MY",ylab="Correlation^2")
  text(years,tbl$iso_cor,rownames(tbl),pos=4)
  abline(lsfit(years,tbl$iso_cor))

  points(years,tbl$acc_cor)
  text(years,tbl$acc_cor,rownames(tbl),pos=1)
  abline(lsfit(years,tbl$acc_cor),lty=2)
  dev.off()
}
