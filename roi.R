writeSubsetROI <- function(set, win.size, fname) {
  old.chrs <- set[,'chr']
  chrs <- unlist(sapply(old.chrs, numToChr))
  out <- data.frame(chr=chrs, start=set[,'pos'], end=set[,'pos'] + win.size)
  write.table(out, file=fname, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
}

loadROI <- function(fname) {
  roi <- read.table(fname, sep='\t', header=FALSE, colClasses=c('character', 'integer', 'integer'))
  roi[,1] <- sapply(roi[,1], chrToNum)
  return(as.matrix(roi))
}

mergeROIWindows <- function(roi) {
  prev.chr <- roi[1, 1]
  prev.start <- roi[1, 2]
  prev.end <- roi[1, 3]
  nrows <- length(roi[,1])
  mtx <- foreach(i = 2:nrows, .combine = rbind) %do% {
    curr.roi <- roi[i,]
    if ((curr.roi[1] == prev.chr) && (curr.roi[2] == prev.end)) {
      prev.end <- curr.roi[3]
      if (i < nrows) {
        return(NA)
      } else {
        return(c(prev.chr, prev.start, prev.end))
      }
    } else {
      to.ret <- c(prev.chr, prev.start, prev.end)
      prev.chr <- curr.roi[1]
      prev.start <- curr.roi[2]
      prev.end <- curr.roi[3]
      return(to.ret)
    }
  }
  colnames(mtx) <- c("chr", "start", "end")
  rownames(mtx) <- 1:length(mtx[,1])
  return(mtx[which(!is.na(mtx[,1])),])
}
