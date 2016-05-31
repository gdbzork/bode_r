eg2bed <- function(rows) {
  require('org.Hs.eg.db')
  left <- as.data.frame(org.Hs.egCHRLOC)[rows,]
  right <- as.data.frame(org.Hs.egCHRLOCEND)[rows,]
  tss <- left$start_location
  tss[tss < 0] <- abs(right$end_location[tss < 0])
  strand <- vector(mode="character",length=dim(left)[1])
  strand[left$start_location > 0] <- '+'
  strand[left$start_location < 0] <- '-'
  chrom <- paste("chr",left$Chromosome,sep="")
  name <- paste(chrom,"_",tss,"_",left$gene_id,sep="")
  bed <- data.frame(chrom,tss,tss+1,name,0,strand)
  colnames(bed) <- c("chrom","left","right","name","score","strand")
  return(bed)
}
