overlapOne = function(chr,st,en,target) {
  filt <- chr = target$chr &&
          && ((st >= target$start && st <= target$end)
              || (en >= target$start && en <= target$end))
  return sum(filt) > 0
}

isContained <- function(thisLeft,thisRight,thatLeft,thatRight) {
  startMatch <- (thisLeft >= thatLeft) & (thisLeft <= thatRight)
  endMatch <- (thisRight >= thatLeft) & (thisRight <= thatRight)
  return (startMatch | endMatch)
}

overlapOne <- function(src,target) {
  chrMatch <- src[1] == target$chr
  srcLeft = as.numeric(src[2])
  srcRight = as.numeric(src[3])
  startMatch <- isContained(srcLeft,srcRight,target$start,target$end)
  endMatch <- isContained(target$start,target$end,srcLeft,srcRight)
  filt <- chrMatch & (startMatch | endMatch)
  return(sum(filt))
}

overlap = function(src,target) {
  filt <- apply(src,1,overlapOne,target)
  return (filt)
}
