# count up distinct anticodons from a list of tRNA names
spfunc = function(x) { strsplit(x,"_",fixed=TRUE)[[1]][2] }
anticodons = function(vec) { 
  x = apply(as.matrix(vec),1,spfunc)
  levels(as.factor(x))
}

# Plots limits of agreement as described by Bland & Altman (Statistical Methods in Medical Research vol 8 no. 2 p. 135)
# x is the vector of average measurements, y the vector of differences
# limix, limiy are the x and y plotting limits
slopelimits <- function(x,y,limix,limiy){
#   limix <- c(min(x),max(x))
#   limiy <- c(min(y),max(y))
    model <- lm(y~x)
    # test to see whether the limits may be considered as horizontal lines
    if(summary(model)$coefficients[2,4] < 0.05){
        model2 <- lm(abs(residuals(model))~x)
        # test to see whether the size of the error depends on the true value
        if(summary(model2)$coefficients[2,4] < 0.05){
            line.central <- c(limix[1],model$coefficients[1]+ limix[1]*model$coefficients[2], limix[2], model$coefficients[1] + limix[2]*model$coefficients[2])
            line.lower <- c(limix[1],model$coefficients[1] +limix[1]*model$coefficients[2] -2.46*(model2$coefficients[1]+limix[1]*model2$coefficients[2]), limix[2], model$coefficients[1] + limix[2]*model$coefficients[2]-2.46*(model2$coefficients[1]+limix[2]*model2$coefficients[2]))
            line.upper <- c(limix[1],model$coefficients[1] +limix[1]*model$coefficients[2] +2.46*(model2$coefficients[1]+limix[1]*model2$coefficients[2]), limix[2], model$coefficients[1] + limix[2]*model$coefficients[2]+2.46*(model2$coefficients[1]+limix[2]*model2$coefficients[2]))
                }
        else{
            line.central <- c(limix[1],model$coefficients[1] +limix[1]*model$coefficients[2], limix[2], model$coefficients[1] + limix[2]*model$coefficients[2])
            line.lower <- c(limix[1],model$coefficients[1] +limix[1]*model$coefficients[2] -1.96*summary(model)$sigma, limix[2], model$coefficients[1] + limix[2]*model$coefficients[2]-1.96*summary(model)$sigma)
            line.upper <- c(limix[1],model$coefficients[1] +limix[1]*model$coefficients[2] +1.96*summary(model)$sigma, limix[2], model$coefficients[1] + limix[2]*model$coefficients[2]+1.96*summary(model)$sigma)  
            print(cat("95% reference interval is","\n"))
            print(cat("lower limit:", round(model$coefficients[1] - 1.96*summary(model)$sigma,digits=3),"+ (",round(model$coefficients[2],digits=3),"* a)",fill=T))     
            print(cat("upper limit:", round(model$coefficients[1] + 1.96*summary(model)$sigma,digits=3),"+ (",round(model$coefficients[2],digits=3),"* a)","\n","where a is the magnitude of the average","\n", fill=F))
    }   
}
    else{
        line.central <- c(limix[1],mean(y), limix[2], mean(y))
        line.lower <- c(limix[1],mean(y) -1.96*sqrt(var(y)), limix[2], mean(y)-1.96*sqrt(var(y)))
        line.upper <- c(limix[1],mean(y) +1.96*sqrt(var(y)), limix[2], mean(y)+1.96*sqrt(var(y)))
        se.lim <- ((1.71)/sqrt(length(x)))*sd(y)
        print("95% reference interval is")
        print(c(round(line.lower[2],digits=3),round(line.upper[2],digits=3)))
        print("95% CI for reference limits are")
        print(cat("lower limit: ",round(line.lower[2]-1.96*se.lim,digits=3),"to",round(line.lower[2]+1.96*se.lim,digits=3),fill=T))
        print(cat("upper limit: ",round(line.upper[2]-1.96*se.lim,digits=3),round(line.upper[2]+1.96*se.lim,digits=3),fill=T))
    }
    # plot the graph
    limiy <- c(min(limiy[1],line.lower[c(2,4)]),max(limiy[2],line.upper[c(2,4)]))
    smoothScatter(x,y,ylim=limiy,xlim=limix)
    segments(line.central[1],line.central[2],line.central[3],line.central[4])
    segments(line.lower[1],line.lower[2],line.lower[3],line.lower[4],lty=2)
    segments(line.upper[1],line.upper[2],line.upper[3],line.upper[4],lty=2)
}

# smoothPairs -- pairs plot with scatterSmooth and "limits of agreement" plots
# "copy" is a matrix of values to be plotted, column against column.
smoothPairs = function(copy) {
  par(mfrow=c(ncol(copy),ncol(copy)),oma=c(0,0,5,0),mar=c(1,1,1,1))
  for (a in 1:ncol(copy)) {
    for (b in 1:ncol(copy)) {
      if (a != b) {
        if (a < b) {
          smoothScatter(copy[,b],copy[,a],xlab=colnames(copy)[b],ylab=colnames(copy)[a],nrpoints=100)
        } else {
           y = (copy[,a] + copy[,b]) / 2
           z = copy[,a] - copy[,b]
           slopelimits(y,z,c(0,10),c(-10,10))
        }
      } else { 
        plot(1, type="n", axes=F, xlab="", ylab="")
        text(1,1,colnames(copy)[a],cex=1.5)
      }
    }
  }
}

# odds and ends for array significance.
upreg <- function(val,mean,stdev,thresh) { val > mean + stdev * thresh }
downreg <- function(val,mean,stdev,thresh) { val < mean - stdev * thresh }
outreg <- function(val,mean,stdev,thresh) { abs(val - mean) > stdev * thresh }

outliers <- function(values,threshold) {
  v.sd <- apply(values,1,sd)
  v.mean <- apply(values,1,mean)
  v.outlier <- outreg(values,v.mean,v.sd,threshold)
  v.outsum <- apply(v.outlier,1,sum)
  values[v.outsum > 0,]
}

makeRaggedList <- function(data,labels) {
  l <- list()
  for (lab in labels) {
    x <- data[data[[1]] == lab,]
    y <- list(x[[2]])
    names(y) <- c(lab)
    l <- c(l,y)
  }
  l
}

makeBoxplot <- function(data,labels,title) {
  boxplot(data,range=1.0,names=labels,main=title,ylim=c(0,2300))
  ms <- sapply(data,mean)
  points(ms)
  ma <- approx(ms,xout=c(1,2,3,4))
  lines(ma)
}

makeMultiLineChart <- function(data,ylim,xlab,ylab,main) {
  pcount = dim(data)[1]
  lcount = dim(data)[2]
  xbar = 1:pcount
  xb = xbar[!is.na(data[,1])]
  db = data[!is.na(data[,1]),1]
  plot(xb,db,ylim=c(0,ylim),type='l',col=1,lwd=2,xaxt='n',xlab=xlab,ylab=ylab,main=main)
  points(xb,db,col=1,pch=1,lwd=2)
  for (i in 2:lcount) {
    xb = xbar[!is.na(data[,i])]
    db = data[!is.na(data[,i]),i]
    lines(xb,db,ylim=c(0,ylim),col=i,lwd=2)
    points(xb,db,col=i,pch=i,lwd=2)
  }
  axis(1,at=xbar,labels=rownames(data))
  legend('bottomright',legend=colnames(data),inset=0.05,lwd=2,col=1:lcount,pch=1:lcount)
}
    
