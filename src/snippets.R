library(IRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(DiffBind)

macs2granges <- function(filename) {
  df <- read.table(filename,blank.lines.skip=TRUE,
                            header=TRUE,
                            comment.char="#",
                            sep="\t",
                            stringsAsFactors=FALSE)
  gr <- GRanges(seqnames=Rle(df$chr),
                ranges=IRanges(df$start,df$end),
                strand=Rle("*",nrow(df)),
                score=df$fold_enrichment)
  return(gr)
}

dba2granges <- function(d) {
  gr <- GRanges(seqnames=Rle(d$chrmap[d$vectors$CHR]),
                ranges=IRanges(d$vectors$START,d$vectors$END),
                strand=Rle("*",nrow(d$vectors)))
  return(gr)
}

