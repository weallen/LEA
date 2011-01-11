dsetToPath <- function(dsetname) {
  return(paste("data", paste(dsetname, "bed", sep="."), sep="/"))
}

rescale.vector <- function(chip.data, scale.factor) {
  res.len <- ceiling(length(chip.data) / scale.factor)
  in.len <- length(chip.data)
  rescaled <- matrix(0, res.len)
  for (i in (1:res.len)) {
    start.idx <- (i-1)*scale.factor
    end.idx <- min(i*scale.factor, in.len)
    rescaled[i] <- sum(chip.data[start.idx:end.idx])
  }
  return(rescaled)
}

loadROIFile <- function(roi.file.name) {
  if (file.exists(roi.file.name)) {
    roi <- read.table(roi.file.name, header=F)
  } else {
    stop("Can't find ROI file")
  }
  colnames(roi) <- c("chr", "start", "stop", "name")
  if (length(unique(roi$chr)) == 1) {
    roi.chr <- unique(roi$chr)
  } else {
    roi.chr <- mixedsort(unique(roi$chr))
  }
  return(list(roi=roi, chr=roi.chr))
}

numToChr <- function(chr) {
  if (chr == 19.0) 
    return("chrX")
  if (chr == 20.0)
    return("chrY")
  if (chr == 21.0)
    return("chrM")
  else 
    return(paste("chr", chr, sep=""))
}

