fn2species <- function(fn) {
  strsplit(fn,"_")[[1]][1]
}

fn2tissue <- function(fn) {
  strsplit(fn,"_")[[1]][2]
}

mergeCounts <- function(fns,keyfn) {
  tags = sapply(fns,keyfn)
  data = c()
  for (fn in fns) {
    tb = read.table(fn,header=FALSE)
    data = cbind(data,tb$V1)
  }
  colnames(data) = tags
  rownames(data) = tb$V2
  data
}

ourHeatMap = function(fns,keyfn,main,ylab,xlab) {
  library(preprocessCore)
  data = mergeCounts(fns,keyfn)
  datn = normalize.quantiles(data)
  rownames(datn) = rownames(data)
  colnames(datn) = colnames(data)
  sq = seq(from=199,to=0)
  heatmap(datn,ylab=ylab,xlab=xlab,scale="column",col=rgb(1,sq/199,sq/199),main=main)
  datn
}

multiPlot = function(tbl) {
  rco = dim(tbl)[2]
  cnames = colnames(tbl)
  rnames = rownames(tbl)
  par(mfrow=c(rco,rco),mar=c(1,1,1,1))
  sq = 1:rco
  mval = max(tbl)
  for (i in sq) {
    for (j in sq) {
      if (i==j) {
        plot(0, axes = TRUE, type = "n", tcl = -0, col.axis = "white", ylab = "", xlab= "")
        text(1,0,cnames[i],pos=3,cex=1.5)
      } else if (i<j) {
        plot(tbl[,j],tbl[,i],xlim=c(0,mval),ylim=c(0,mval),ylab="",xlab="")
        text(tbl[,j],tbl[,i],rnames,pos=4,offset=0.5)
        abline(lsfit(tbl[,j],tbl[,i]))
      } else {
        plot(0, axes = TRUE, type = "n", tcl = -0, col.axis = "white", ylab = "", xlab= "")
        dcor = cor(tbl[,j],tbl[,i],method="spearman")
        dcorsq = dcor^2
        text(1,0,dcorsq,pos=3)
      }
    }
  }
}
