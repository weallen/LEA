batchMode = TRUE
if (batchMode) {
  cmdArgs = commandArgs(TRUE)
  print(cmdArgs)
}


GENOME.STATS <- "/gpfs/home/wallen/data/genome_data/mm9_window_1kb_stats.txt"

compareCpGAndCpAMeth <- function() {
#  setwd("/gpfs/home/wallen/data/genome_data")
#  window.1kb.stats <- read.table("mm9_window_1kb_stats.txt",
#                                 header=TRUE, sep="\t")
  setwd("~/experiment/experiment/stavros_data")
  load("/gpfs/home/wallen/data/genome_data/mm9_window_1kb_stats.Rdata")
  omp.medip <- load.DipData("omp_medip_1kb")
  omp.hmedip <- load.DipData("omp_hmedip_1kb")
  cat("Sampling...\n")
  idx <- sample(1:length(omp.medip$genome.data[,'raw']), 1000000)
  png("cpn_corr.png")
  par(mfrow=c(2,2))
  cat("Drawing plot...\n")
  plot(mm9.window.1kb.stats$CpA[idx],omp.medip$genome.data[idx, 'raw'],
       pch=20, cex=.1, col=rgb(0, 0, 100, 50, maxColorValue=255),
       ylim=c(0, 100), main="CpA meC")
  plot(mm9.window.1kb.stats$CpG[idx],omp.medip$genome.data[idx, 'raw'],
       pch=20, cex=.1, col=rgb(0, 0, 100, 50, maxColorValue=255),
       ylim=c(0, 100), main="CpG meC")
  plot(mm9.window.1kb.stats$CpA[idx],omp.hmedip$genome.data[idx, 'raw'],
       pch=20, cex=.1, col=rgb(0,0,100, 50, maxColorValue=255),
       ylim=c(0, 100), main="CpA hmeC")
  plot(mm9.window.1kb.stats$CpG[idx],omp.hmedip$genome.data[idx, 'raw'],
       pch=20, cex=.1, col=rgb(0,0,100,50, maxColorValue=255),
       ylim=c(0, 100), main="CpG hmeC")
  dev.off()
}

load("/gpfs/home/wallen/data/genome_data/mm9_window_1kb_stats.Rdata")
omp.medip <- load.DipData("omp_medip_1kb")$genome.data[,'raw']
  
