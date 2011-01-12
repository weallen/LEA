library(bigmemory)
library(multicore)
library(doMC)

options(bigmemory.typecast, warning=FALSE)
options(bigmemory.allow.dimnames=TRUE)
registerDoMC(cores=2)


DipData <- function(data.set.name, genome.chr=NULL, genome.pos=NULL,
                    genome.raw=NULL, genome.norm=NULL, chr.names=NULL,
                    chr.lengths=NULL, bin.size=NULL) {
  obj <- list(name = data.set.name,
              genome.chr = genome.chr,
              genome.pos = genome.pos,
              genome.raw = genome.raw,
              genome.norm = genome.norm,
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
  .writeDipDataCol(dip.data$genome.chr, data.path, "genome_chr")
  .writeDipDataCol(dip.data$genome.pos, data.path, "genome_pos")
  .writeDipDataCol(dip.data$genome.raw, data.path, "genome_raw")
  .writeDipDataCol(dip.data$genome.norm, data.path, "genome_norm")
}

.writeDipDataCol <- function(dip.data.col, data.path, colname) {
  write.big.matrix(dip.data.col, paste(data.path, paste(colname, "txt", sep="."), sep="/"))
}

.readDipDataCol <- function(data.path, colname, dtype) {
  col.path <- paste(colname, "txt", sep=".")
  col.bin.path <- paste(colname, "bin", sep=".")
  return(read.big.matrix(paste(data.path, col.path, sep="/"),
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
  
  obj$genome.chr <- .readDipDataCol(data.path, "genome_chr", "char")
  obj$genome.pos <- .readDipDataCol(data.path, "genome_pos", "double")
  obj$genome.raw <- .readDipDataCol(data.path, "genome_raw", "double")
  obj$genome.norm <- .readDipDataCol(data.path, "genome_norm", "double")
  chr.data <- read.table(paste(data.path, dip.data.name, sep="/"), header=TRUE)
  obj$chr.names <- chr.data$chr.names
  obj$chr.lengths <- chr.data$chr.lengths
  obj$bin.size <- chr.data$bin.size
  return(obj)
}

diffEnrichment.DipData <- function(dd1, dd2) {
  apply
}
