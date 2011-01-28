# TODO Fix merging
source("util.R")

writeSubsetROI <- function(set, win.size, fname, merge=FALSE) {
  old.chrs <- set[,'chr']
  chrs <- unlist(sapply(old.chrs, numToChr))
  out <- data.frame(chr=chrs, start=set[,'pos'], end=set[,'pos'] + win.size)
  if (merge) {
    out <- mergeROIWindows(out)
  }
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

diffEnrichmentROI <- function(enrich1, enrich2) {
  roi1.sum <- sum(enrich1)
  roi2.sum <- sum(enrich2)
  p.val <- foreach(i=icount(length(enrich1)), .combine=c) %do% {    
    cont.table <- matrix(c(enrich1[i], roi1.sum - enrich1[i],
                           enrich2[i], roi2.sum - enrich2[i]), nrow=2, ncol=2)
    return(fast.fisher(cont.table)$p.val)
  }
  return(p.val)
}

saveSignificantROI <- function(in.roi, p.vals, path, cutoff=0.1) {
  in.roi <- as.data.frame(in.roi)
  cat("Found", sum(p.vals < cutoff), "roi less than threshold", cutoff, "\n")
  if (sum(p.vals < cutoff) > 0) {
    subset.roi <- in.roi[which(p.vals < cutoff),]
    chrs <- sapply(subset.roi$chr, numToChr)
    subset.roi$chr <- chrs
    write.table(subset.roi, file=path, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  }
}
