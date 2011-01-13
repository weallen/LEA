library(MEDIPS)
library(gtools)
library(edgeR)
library(plyr)
library(multtest)
library(BSgenome.Mmusculus.UCSC.mm9)
library(GDD)
library(bigmemory)

setwd("/gpfs/home/wallen/experiment/experiment/stavros_data")

source("~/src/LEA/medips.R")
source("~/src/LEA/util.R")

convertMedipsToDipData <- function() {
  cat("omp\n")
  load("omp_hmedip.Rdata")
  saveMedipsForDipData(omp.hmedip, "omp_hmedip")
  rm(omp.hmedip)

  cat("ngn\n")
  load("ngn_hmedip.Rdata")
  saveMedipsForDipData(ngn.hmedip, "ngn_hmedip")
  rm(ngn.hmedip) 

  cat('icam\n')
  load("icam_hmedip.Rdata")
  saveMedipsForDipData(icam.hmedip, "icam_hmedip")
  rm(icam.hmedip)

  cat('moe\n')
  load("moe_hmedip.Rdata")
  saveMedipsForDipData(moe.hmedip, "moe_hmedip")
  rm(moe.hmedip)

  cat('moe_ac3\n')
  load("moe_ac3_hmedip.Rdata")
  saveMedipsForDipData(moe.ac3.hmedip, "moe_ac3_hmedip")
  rm(moe.ac3.hmedip)

  gc()
}

rescaleAndDiff <- function(medips1, medips2, window.size=1000) {
  cat("Rescaling medips1\n")
  medips1.res <- reduceMedipsResolution(medips1, window.size)
  cat("Rescaling medips2\n")
  medips2.res <- reduceMedipsResolution(medips2, window.size)
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

selectSignificants <- function(medips.diff, fdr.cutoff=0.01) {
  library(multtest)
  q.values <- mt.rawp2adjp(medips.diff$p.value.wilcox, proc="BH")
  
}
# SCRIPT

files <- c("omp_hmedip.bed", "ngn_hmedip.bed", "icam_hmedip.bed", "moe_hmedip.bed", "moe_ac3_hmedip.bed")
chrs <- c(paste("chr", 1:19, sep=""), "chrX", "chrY")
# HMEDIP
#loadAndSaveAllHmeDipData()

#moe.ac3.diff <- MEDIPS.methylProfiling(data1=moe.hmedip, data2=moe.ac3.hmedip, select=1, frame_size=10000)
#moe.ac3.sig <- MEDIPS.mergeFrames(frames=MEDIPS.selectSignificants(frames=moe.ac3.diff, control=F, p.value=1e-04))
