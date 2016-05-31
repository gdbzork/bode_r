#!/usr/bin/Rscript

percentile = function(x,prob=0.75) {
  s = sort(x)
  pos = round(length(s) * prob)
  s[pos]
}

plotDensities = function(configfile,title,show_var=FALSE,func=mean,vfunc=sd) {
  config = read.table(configfile,sep='\t',as.is=TRUE)
  colnames(config) = c('color','name','readcount','datafile')

  meanlist <- list()
  sdlist = list()
  maxvalues <- vector()
  xaxis_max = 0
  maxcount = max(config$readcount)
  
  for (i in seq(from=1,to=dim(config)[1])) {
    zork = config$datafile[i]
    dtable = read.table(zork)
    scale_factor = maxcount / config$readcount[i]
#    dtable = dtable * scale_factor
    dmeans = apply(dtable,2,func)
    maxvalues = c(maxvalues,max(dmeans))
    meanlist = c(meanlist,list(dmeans))
    if (show_var) {
      dsds = apply(dtable,2,vfunc)
      maxvalues = c(maxvalues,max(dsds))
      sdlist = c(sdlist,list(dsds))
    }
    print(c(zork,' scale=',scale_factor),quote=FALSE)
  }

  xlab = "Base position relative to 5' end of tRNA"
  ylab = "Median Read Coverage"
  realmax = max(maxvalues)

  xaxis_max = length(meanlist[[1]])
  xaxis_middle = xaxis_max / 2
  xaxis = seq(from=-xaxis_middle,to=xaxis_middle-1)
  plot(xaxis,meanlist[[1]],col=config$color[1],type='l',ylim=c(0,realmax),xlab=xlab,ylab=ylab,main=title)
  if (show_var) {
    points(xaxis,sdlist[[1]],type='l',lty=3,col=config$color[1])
  }
  for (i in seq(from=2,to=dim(config)[1])) {
    points(xaxis,meanlist[[i]],col=config$color[i],type='l')
    if (show_var) {
      points(xaxis,sdlist[[i]],type='l',lty=3,col=config$color[i])
    }
  }

  legend('topright',legend=config$name,col=config$color,lwd=2)
}
