library(bigmemory)
library(bigtabulate)
library(biganalytics)
library(bigalgebra)
library(bbmle)
library(MASS)
library(multicore)
library(foreach)
library(doMC)
#library(qvalue)

#options(bigmemory.typecast, warning=FALSE)
options(bigmemory.allow.dimnames=TRUE)
registerDoMC(cores=10)

source("~/src/LEA/util.R")

DipData <- function(data.set.name, genome.data=NULL, chr.names=NULL,
                    chr.lengths=NULL, bin.size=NULL) {
  # NOTE: columns in genome.data are c('chr', 'pos', 'raw', 'norm')
  obj <- list(name = data.set.name,
              genome.data = genome.data,
              chr.names = chr.names,
              chr.lengths = chr.lengths,
              bin.size = bin.size)                
  class(obj) <- "DipData" 
  return(obj)
}

normDist.DipData <- function(dd1, dd2) {
  dd1.data <- as.vector(dd1$genome.data[,'norm'])
  dd2.data <- as.vector(dd2$genome.data[,'norm'])
  return(sqrt(sum(dd1.data - dd2.data)**2))
}

getNormEnrich.DipData <- function(dd) {
  return(dd$genome.data[,'norm'])
}

save.DipData <- function(dip.data) {
  data.path <- paste(dip.data$name, "dipdata", sep=".")
  if (file.exists(data.path)) {
    unlink(data.path, recursive=TRUE)
  }
  dir.create(data.path)
  write.table(dip.data[c('name', 'chr.names', 'chr.lengths')], file=paste(data.path, dip.data$name, sep="/"), sep=",")
  .writeDipData(dip.data$genome.chr, data.path, "genome_data")
}

.writeDipData <- function(dip.data, data.path, colname) {
  write.big.matrix(dip.data.col, paste(data.path, paste(colname, "txt", sep="."), sep="/"))
}

.readDipData <- function(colname, dtype) {
  col.path <- paste(colname, "txt", sep=".")
  col.bin.path <- paste(colname, "bin", sep=".")
  return(read.big.matrix(col.path, header=TRUE, 
                         type=dtype, backingfile=col.bin.path,
                         backingpath=".",
                         descriptorfile=paste(colname, "desc",sep=".")))
}

load.DipData <- function(dip.data.name) {
  .doLoad <- function() {
    obj <- DipData(dip.data.name)
    if (!file.exists("genome_data.bin") ||
        !file.exists("genome_data.desc")) {
      obj$genome.data <- .readDipData("genome_data", "double")    
    }
    else {
      obj$genome.data <- attach.big.matrix(dget("genome_data.desc"))
    }
    chr.data <- read.table(dip.data.name, header=TRUE, sep=",")
    obj$chr.names <- chr.data$chr.names
    obj$chr.lengths <- chr.data$chr.lengths
    obj$bin.size <- chr.data$bin.size
    return(obj)
  }
  # kind of a kludge but w/e
  old.wd <- getwd()
  data.path <- paste("dipdata",paste(dip.data.name, "dipdata", sep="."), sep="/")
  if (!file.exists(data.path)) {
    stop("Dipdata file does not exist")
  }
  setwd(data.path) 
  tryCatch(obj <- .doLoad(), finally=setwd(old.wd))
  return(obj)
}

foldChange.DipData <- function(dd1, dd2, type="raw") {
  require(gtools)
  if (type != "raw" && type != "norm") {
    stop("Type must be raw or norm")
  }
  chr.idx <- bigsplit(dd$genome.data, 'chr')
  fc <- foreach(chr=chr.idx, .combine=c) %dopar% {
    sapply(chr, function(i) {
      return(foldchange(dd1$genome.data[i, type], dd2$genome.data[i, type]))
    })
  }
  return(fc)
}

diffEnrichment.DipData <- function(dd1, dd2) {
  chr.idx <- bigsplit(dd1$genome.data, 'chr')
  freq <- list(dd1=dd1$genome.data, dd2=dd2$genome.data)
  dd1.sum <- colsum(dd1$genome.data, "raw")
  dd2.sum <- colsum(dd2$genome.data, "raw")
  p.val <- foreach(chr=chr.idx, .combine=c) %dopar% {    
    sapply(chr,  function(i) {
      cont.table <- matrix(c(dd1$genome.data[i, "raw"], dd1.sum - dd1$genome.data[i, "raw"],
                             dd2$genome.data[i, "raw"], dd2.sum - dd2$genome.data[i, "raw"]), nrow=2, ncol=2)
      return(fast.fisher(cont.table)$p.val)
    })
  }
  return(p.val)
}

# ROI is data.frame(chr=, start=, end=)
subsetByROI.DipData <- function(dd, in.roi) {
  bin.size <- dd$bin.size[1]
  cat(bin.size, '\n')
  absolute.pos <- c(1, cumsum(ceiling(dd$chr.lengths / dd$bin.size[[1]])))
  idx <- foreach(curr.roi = isplitRows(in.roi, chunkSize=1), .combine=c) %do% {
    curr.chr <- curr.roi[[1]]
    curr.chr.start <- absolute.pos[curr.chr]
    curr.roi.start <- curr.chr.start + ceiling(curr.roi[[2]] / bin.size)
    cat(curr.roi.start,'\n')
    curr.roi.end <- curr.chr.start + ceiling(curr.roi[[3]] / bin.size)
    cat(curr.roi.end,'\n')
    return(seq(curr.roi.start,curr.roi.end))
  }
  return(DipData(dd$name, dd$genome.data[idx,], dd$chr.names, dd$chr.lengths, dd$bin.size))
}

computeAndSaveDiffEnrich.DipData <- function(dd1, dd2, fname, roi=NULL) {
  cat("calling diff enrich", fname, '\n')
  p.val <- diffEnrichment.DipData(dd1, dd2)
  cat("adjusting  pvals", fname, '\n')
  q.val <- callSignificant(p.val)
  tab <- data.frame(chr=dd1$genome.data[,1], pos=dd1$genome.data[,2], qval=q.val)
  cat('writing data', fname, '\n')
  write.table(tab, file=paste("diff_enrich", fname, sep="/"), quote=FALSE, sep=",", row.names=FALSE)
}

# Coverts matrix of differentially enriched windows to matrix of spans
# of differential enrichment (merging nearby windows)
# uses a greedy algorithm for merging, allowing for missing windows up to some
# cutoff
mergeDiffEnrich.DipData <- function(dd, p.vals, p.cutoff=0.05, gaps.cutoff=10) {
  diff.enriched <- dd$genome.data[p.vals < p.cutoff]
  bin.size <- dd$bin.size[1]
  chr.idx <- bigsplit(diff.enriched, "chr")
  mtx <- matrix(0, length(diff.enriched[,"chr"]), 3)
  mtx.idx <- 0
  foreach(idx=icount(length(diff.enriched[,"chr"]))) %do% {
    curr.chr <- diff.enriched[idx, "chr"]
    curr.span.start <- diff.enriched[idx, "pos"]
    curr.span.end <- curr.span.start + bin.size
    pos <- diff.enriched[idx, "pos"]
    if (pos == curr.span.end) {
      curr.span.end <- pos + bin.size
    }
    else if (pos > curr.span.end && (pos <= (curr.span.end + (bin.size * (gaps.cutoff + 1))))) {
      curr.span.end <- pos + bin.size
    }
    else {
      curr.span.start <- pos
      curr.span.end <- pos + bin.size
      mtx[mtx.idx, 1] <- curr.chr
      mtx[mtx.idx, 2] <- curr.span.start
      mtx[mtx.idx, 3] <- curr.span.end
      mtx.idx <- mtx.idx + 1
    }
  }
  return(mtx[mtx[,1] > 0,])
}


# Usually you'd want to do multiple hypothesis adjustment before doing this
saveDiffEnrichedWIG.DipData <- function(dd, fname, p.vals, cutoff=0.05) {
  diff.enriched <- dd$genome.data[p.vals < cutoff]
  chr.idx <- bigsplit(diff.enriched, "chr")
  write("track type=wiggle_0", file=fname)
  for(chr in dd$chr.names) {
    curr.enriched <- diff.enriched[which(diff.enriched[,"chr"] == chrToNum(chr)), ]
    if (length(curr.enriched) > 0) {
      chr.start <- min(curr.enriched[,"pos"])
      chrom.str <- paste("fixedStep chrom=",chr," start=", chr.start, " step=100 span=100", sep="") 
      write(chrom.str, file=fname, append=TRUE)
      write(curr.enriched[,"pos"], file=fname, append=TRUE, sep="\n")
    }
  }
}

getDiffEnrichPos.DipData <- function(dd, qvals, fdr=0.05) {
  return(subset.DipData(dd, which(qvals < .05)))
}

subset.DipData <- function(dd, positions) {  
  return(dd$genome.data[positions,])
}

# Fits a normal distribution to the log2 normalized data
thresholdParams <- function(data) {
  fixed.data <- data[data > -Inf]
  fit <- fitdistr(fixed.data, 'normal')
  return(list(mean=fit$estimate['mean'], sd=fit$estimate['sd']))
}

# XXX not yet fully implemented
normThreshold.DipData <- function(dd, fdr = 0.1) {
  data <- log2(dd$genome.data[dd$genome.data[,'norm']>0,'norm'])
  params <- thresholdParams(data)
  z.vals <- (data - params$mean) / params$sd

  ### NOTE RETURNS CONSTANT
  return(0.5)
}

callEnriched.DipData <- function(dd) {
  thresh <- normThreshold.DipData(dd)
  return(dd$genome.data[dd$genome.data[,'norm'] > thresh, ])
}


writeSubsetROI <- function(set, win.size, fname) {
  chrs <- sapply(set[,'chr'], numToChr)
  out <- data.frame(chr=numToChr(set[,'chr'], start=set[,'pos'], end=set[,'pos'] + win.size)
  write.table(out, file=fname, sep='\t', header=FALSE, quote=FALSE, row.names=FALSE)
}



