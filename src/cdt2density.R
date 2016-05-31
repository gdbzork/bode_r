loadCDT <- function(fn) {
}

cdt2density <- function(cdt) {
  d <- apply(cdt,c(2),mean)
  q1 <- 
}

fn = commandArgs(T)

files = read.table(fn,header=F)
colnames(files) = c('tag','fn')
