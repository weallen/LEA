library(bigmemory)
library(bigtabulate)
library(biganalytics)
library(multicore)
library(foreach)
library(doMC)

#options(bigmemory.typecast, warning=FALSE)
options(bigmemory.allow.dimnames=TRUE)
registerDoMC(cores=16)


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

.readDipData <- function(data.path, colname, dtype) {
  col.path <- paste(colname, "txt", sep=".")
  col.bin.path <- paste(colname, "bin", sep=".")
  return(read.big.matrix(paste(data.path, col.path, sep="/"), header=TRUE, 
                         type=dtype, backingfile=col.bin.path,
                         backingpath=data.path,
                         descriptorfile=paste(colname, "desc",sep=".")))
}

load.DipData <- function(dip.data.name) {
  data.path <- paste(dip.data.name, "dipdata", sep=".")
  if (!file.exists(data.path)) {
    stop("Error: Dipdata file does not exist")
  }
  obj <- DipData(dip.data.name)
  obj$genome.data <- .readDipData(data.path, "genome_data", "double")
  chr.data <- read.table(paste(data.path, dip.data.name, sep="/"), header=TRUE)
  obj$chr.names <- chr.data$chr.names
  obj$chr.lengths <- chr.data$chr.lengths
  obj$bin.size <- chr.data$bin.size
  return(obj)
}

foldChange.DipData <- function(dd1, dd2) {
  
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
  
#fisherTestResults = apply(frequencyTable,1,function(X)
#  { contingencyTable = matrix(c(X["sampleReads"],sampleTotal-X["sampleReads"],
#                                X["otherSampleReads"],otherSampleTotal-X["otherSampleReads"]),nrow
#                                                                   =2,ncol=2); fast.fisher(contingencyTable) })
#pValueVsOthers = pvalues = unlist(lapply(fisherTestResults, function(X) { X[["p.value"]] }))
#qvalueObject = try(qvalue(signif(pvalues,8))) # "qvalue fails if values are marginally higher than 1, also revert to Bonferroni if FDR estimation fails due to low sample size and odd p-value distri
  
#  fisher.fast
}
