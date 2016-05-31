library(RColorBrewer)

makePlot <- function(means,caption,xlab,ylab,names) {
  lm <- length(means[[1]])
  colpat <- brewer.pal(length(means),'Set1')
  mtop <- max(sapply(means,max))
  coords <- (-lm/2) : (lm/2-1)
  plot(coords,means[[1]],type='l',col=colpat[1],ylim=c(0,mtop+10),xlab=xlab,ylab=ylab)
  for (i in 2:length(means)) {
    lines(coords,means[[i]],type='l',col=colpat[i])
  }
  title(main=caption)
  legend('topright','Legend',names,col=colpat,inset=0.05,lwd=2)
}

mfunc <- function(mat) {
  apply(mat,c(2),mean)
}

makeMeans <- function(raw) {
  means <- lapply(raw,function(mat) { apply(mat,c(2),mean) })
  means
}

processArgs <- function() {
  fns = list()
  for (a in commandArgs(trailingOnly=TRUE)) {
    pair <- strsplit(a,'=')[[1]]
    if (length(pair) == 1) {
      fns = c(fns,list(pair))
    } else {
      assign(pair[1],pair[2],inherits=TRUE)
    }
  }
  fns
}

cap <- "Default Caption"
xlab <- "Default X-Axis Label"
ylab <- "Default Y-Axis Label"
nameset <- NULL
target <- "densityPlot.pdf"

files <- processArgs()
if (!is.null(nameset)) {
  names = strsplit(nameset,',')[[1]]
} else {
  names = lapply(files,function(x) { substr(x,1,5) })
}

data <- list()
for (fn in files) {
  d <- read.table(fn,header=F,stringsAsFactors=F)
  data <- c(data,list(d))
}

mset <- makeMeans(data)
pdf(file=target)
makePlot(mset,cap,xlab,ylab,names)
dev.off()
