library(MEDIPS)
library(gtools)
library(edgeR)
library(plyr)
library(multtest)
library(BSgenome.Mmusculus.UCSC.mm9)
library(GDD)

setwd("/gpfs/home/wallen/experiment/experiment/stavros_data")

medipsToGRanges <- function(medips.object) {
  chrs <- c(paste("chr", 1:18, sep=""), "chrX", "chrY")
  which.idx <- genome_chr(medips.object) %in% chrs
  genome.chr <- Rle(genome_chr(medips.object)[which.idx])
  genome.rng <- IRanges(start=genome_pos(medips.object)[which.idx], width=bin_size(medips.object))
#  seqlen <- as.vector(chr_lengths(medips.object)[chr_names(medips.object) %in% chrs], mode='integer'
  return(GRanges(seqnames=genome.chr, ranges=genome.rng,
                 reads=genome_raw(medips.object)[which.idx]))
}

loadMedips <- function(fname) {
  d <- MEDIPS.readAlignedSequences(BSgenome = "BSgenome.Mmusculus.UCSC.mm9", file=fname)
  d <- MEDIPS.genomeVector(data=d, bin_size=100, extend=350)
  d <- MEDIPS.getPositions(data=d, pattern="CG")
  d <- MEDIPS.couplingVector(data=d, fragmentLength=700, func="count")
  d <- MEDIPS.calibrationCurve(data=d)
  d <- MEDIPS.normalize(data=d)
  return(d)
}

loadAndSaveAllHmeDipData <- function() {
  omp.hmedip <- loadMedips(dsetToPath("omp_hmedip"))
  save(omp.hmedip,file= "omp_hmedip.Rdata")
  ngn.hmedip <- loadMedips(dsetToPath("ngn_hmedip"))
  save(ngn.hmedip, file="ngn.hmedip.Rdata")
  icam.hmedip <- loadMedips(dsetToPath("icam_hmedip"))
  save(icam.hmedip, file="icam_hmedip.Rdata")
  
  moe.ac3.hmedip <- loadMedips(dsetToPath("moe_ac3_hmedip"))  
  save(moe.ac3.hmedip, file="moe_ac3_hmedip.Rdata")
  moe.hmedip <- loadMedips(dsetToPath("moe_hmedip"))
  save(moe.hmedip, file="moe_hmedip.Rdata")
}

dsetToPath <- function(dsetname) {
  return(paste("data", paste(dsetname, "bed", sep="."), sep="/"))
}
  
plotCalibrateDip <- function(dip.data, fname, chr) {
  GDD(file=fname, type="png", width=1200, height=900)
  MEDIPS.plotCalibrationPlot(data=dip.data, linear=T, plot_chr=chr)
  dev.off()
}

plotCoverageAnalysis <- function(dip.data, fname) {
  GDD(file=fname, type="png", width=1600, height=1200)
  dip.data <- MEDIPS.coverageAnalysis(data=dip.data, extend=350, no_iterations=10)
  MEDIPS.plotCoverage(dip.data)
  dev.off()
  return(dip.data)
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

subsetMedipsByROI <- function(medips.obj, roi) {
  if (!identical(chr_names(medips.obj), roi$chr)) {
    idx <- which(genome_chr(medips.obj) %in% roi$chr)
    if (length(idx) == 0) {
      stop("No chr in common")
    }
    medips.obj@chr_lengths <- chr_lengths(medips.obj)[chr_names(medips.obj) %in% roi$chr]
  } else {
    idx = NULL
  }
  if (!is.null(idx)) {
    medips.obj@genome_chr <- genome_chr(medips.obj)[idx]
    medips.obj@genome_pos <- genome_pos(medips.obj)[idx]
    medips.obj@genome_raw <- genome_raw(medips.obj)[idx]
    medips.obj@genome_norm <- genome_norm(medips.obj)[idx]
    medips.obj@chr_names <- unique(genome_chr(medips.obj))
  }
  return(medips.obj)
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

reduceMedipsResolution <- function(medips.obj, window.size) {

  scale.factor <- window.size / medips.obj@bin_size
  reduced.chr.len <- vector(length=length(medips.obj@chr_names), mode="integer")
  for (i in 1:length(medips.obj@chr_names)) {
    curr.chr <- medips.obj@chr_names[i]
    reduced.chr.len[i] <- ceiling(medips.obj@chr_lengths[i] / window.size)
  }
  scaled.size <- sum(reduced.chr.len)
  chr.pos <- cumsum(reduced.chr.len)
  genome.chr <- vector(scaled.size, mode="integer")
  genome.pos <- vector(scaled.size, mode="integer")
  genome.raw <- vector(scaled.size, mode="integer")
  genome.norm <- vector(scaled.size, mode="numeric")
  medips.obj@bin_size <- window.size
  for (i in 1:length(medips.obj@chr_names)) {
    curr.chr <- medips.obj@chr_names[i]
    if (i == 1) {
      idx.rng <- 1:reduced.chr.len[i]
    } else {
      idx.rng <- (chr.pos[i-1]+1):(chr.pos[i-1]+reduced.chr.len[i])
    }
    genome.chr[idx.rng] <- curr.chr
    genome.pos[idx.rng] <- seq(1, medips.obj@chr_lengths[i], window.size)
    genome.raw[idx.rng] <- rescale.vector(medips.obj@genome_raw[medips.obj@genome_chr == curr.chr], scale.factor)
    genome.norm[idx.rng] <- rescale.vector(medips.obj@genome_norm[medips.obj@genome_chr == curr.chr], scale.factor)
  }
  medips.obj@genome_chr <- genome.chr
  medips.obj@genome_pos <- genome.pos
  medips.obj@genome_raw <- genome.raw
  medips.obj@genome_norm <- genome.norm
  return(medips.obj)
}

# use Fisher's exact test

# fits a poisson GLM to each window TOO SLOW
diffMeth.glm <- function(medips1, medips2, window.size=1000) {
  require(utils)
  num.windows <- length(genome_raw(medips1))
  pb <- txtProgressBar(min=1, max=num.windows, style=3)
  glm.pval <- vector()
  fold.changes <- vector()
  cs <- c(sum(genome_raw(medips1)), sum(genome_raw(medips2)))
  sample.f <- factor(c(1, 2))
  for (i in 1:num.windows) {
    setTxtProgressBar(pb, i)
    s1 <- genome_raw(medips1)[i]
    s2 <- genome_raw(medips2)[i]
    data <- as.vector(unlist(c(s1, s2)))
    glm.curr <- glm(data ~ 1 + sample.f, offset=log(cs), family="poisson")
    glm.pval[i] <- anova(glm.curr, test="Chisq")[5][2,1]
#    fold.changes[i] <- exp(glm.curr$coefficients[1])/(exp(glm.curr$coefficients[1]+glm.curr$coefficients[2]))
  }
  out <- matrix(ncol=2, nrow=num.windows)
  out[,1] <- glm.pval
  out[,2] <- fold.changes
  colnames(out) <- c("pval", "fc")
  return(output)
}


selectSignificants <- function(medips.diff, fdr.cutoff=0.01) {
  library(multtest)
  q.values <- mt.rawp2adjp(medips.diff$p.value.wilcox, proc="BH")
  
}
# SCRIPT

files <- c("omp_hmedip.bed", "ngn_hmedip.bed", "icam_hmedip.bed", "moe_hmedip.bed", "moe_ac3_hmedip.bed")
chrs <- c(paste("chr", 1:18, sep=""), "chrX", "chrY")
# HMEDIP
#loadAndSaveAllHmeDipData()

#moe.ac3.diff <- MEDIPS.methylProfiling(data1=moe.hmedip, data2=moe.ac3.hmedip, select=1, frame_size=10000)
#moe.ac3.sig <- MEDIPS.mergeFrames(frames=MEDIPS.selectSignificants(frames=moe.ac3.diff, control=F, p.value=1e-04))
