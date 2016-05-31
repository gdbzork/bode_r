fn2score <- function(filename,func,tag) {
  t <- read.table(filename)
  d <- apply(t,2,func)
  df <- data.frame(1:length(d),d,tag)
  colnames(df) <- c("base","score","group")
  df$smoothed <- loess(score ~ base, data=df, span=0.005)$fitted
  return(df)
}
  
makeLines <- function(fns,tags) {
  res <- fn2score(fns[1],mean,tags[1])
  for (i in 2:length(fns)) { 
    d <- fn2score(fns[i],mean,tags[i])
    res <- rbind(res,d)
  }
  return(res)
}
