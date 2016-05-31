prop <- function(x) {
  y <- sum(x)
  return (x/y)
}

makeLogo <- function(matrix,outFN="") {
  library(seqLogo)
  if (class(matrix) == "character") {
    df <- read.table(matrix)
  } else {
    df <- matrix
  }
  mat <- as.matrix(df)
  if (dim(mat)[1] != 4) {
    mat <- t(mat)
  }
  if (class(mat[1,1]) == "integer") {
    mat <- apply(mat,2,prop)
  }
  pwm <- makePWM(mat)
  if (outFN != "") {
    pdf(outFN)
  }
  seqLogo(pwm)
  if (outFN != "") {
    dev.off()
  }
}
