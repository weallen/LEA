
library(MEDIPS)
library(gtools)
library(edgeR)
#library(plyr)
library(multtest)
library(BSgenome.Mmusculus.UCSC.mm9)
library(GDD)
library(bigmemory)
library(fields)
library(gplots)

setwd("/gpfs/home/wallen/experiment/experiment/stavros_data")
#setwd("~/Documents/Neuroscience/barnea_lab/rna_seq_experiment/stavros_chip")
source("~/src/LEA/common.R")
source("~/src/LEA/medips.R")
source("~/src/LEA/util.R")
source("~/src/LEA/dipdata.R")

pairwiseDiffEnrich <- function() {
#  pairs <- list(c("moe_hmedip_1kb", "moe_ac3_hmedip_1kb"),
#                c("omp_hmedip_1kb", "icam_hmedip_1kb"),
#                c("moe_medip_1kb", "moe_ac3_medip_1kb"),
#                c("omp_medip_1kb", "icam_medip_1kb"))
  pairs <- list(c("omp_hmedip_1kb", "ngn_hmedip_1kb"),
                c("ngn_hmedip_1kb", "icam_hmedip_1kb"),
                c("omp_medip_1kb", "ngn_medip_1kb"),
                c("ngn_medip_1kb", "icam_medip_1kb"))
  
#  i <- c("moe_medip_1kb", "moe_ac3_medip_1kb")
  foreach(i = pairs) %dopar% {
    fname <- paste(paste(i[1], i[2], sep="_vs_"), "txt", sep=".")
    cat("loading", i[1], '\n')
    dd1 <- load.DipData(i[1])
    cat("loading", i[2], '\n')
    dd2 <- load.DipData(i[2])
    cat('computing and saving diff enrich', fname, '\n')
    computeAndSaveDiffEnrich.DipData(dd1, dd2, fname)
 }                
}

loadMoeMedip <- function() {
  name <- "moe_medip"
  curr.dset.100 <- loadMedips(dsetToPath(name), 100)
  curr.dset.1000 <- loadMedips(dsetToPath(name), 1000)
  curr.dset.10000 <- loadMedips(dsetToPath(name), 10000)
  cat("100bp\n")
  saveMedipsForDipData(curr.dset.100, "moe_medip")
  cat("1kb\n")
  saveMedipsForDipData(curr.dset.1000, "moe_medip_1kb")
  cat("10kb\n")
  saveMedipsForDipData(curr.dset.10000, "moe_medip_10kb")
}

convertAllDipsToDipData25kb <- function() {
  names <- c(MEDIP.DATASETS, HMEDIP.DATASETS)
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading medips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 25000)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, paste(name, "25kb", sep="_"))
  }
}

convertAllDipsToDipData10kb <- function() {
  names <- c("omp_medip", "ngn_medip", "icam_medip", "moe_ac3_medip", "omp_hmedip", "ngn_hmedip", "icam_hmedip", "moe_ac3_hmedip", "moe_hmedip")
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading medips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 10000)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, paste(name, "10kb", sep="_"))
  }
}

convertAllDipsToDipData1kb <- function() {
  names <- c("omp_medip", "ngn_medip", "icam_medip", "moe_ac3_medip", "omp_hmedip", "ngn_hmedip", "icam_hmedip", "moe_ac3_hmedip", "moe_hmedip")
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading medips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 1000)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, paste(name, "1kb", sep="_"))
  }
}

# don't forget to remove random chromosomes before running loadMedips
convertMedipsToDipData <- function() {
  names <- c("omp_medip", "ngn_medip", "icam_medip", "moe_ac3_medip")
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading medips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 100)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, name)
  }
}

convertHMedipsToDipData <- function() {
  names <- c("omp_hmedip", "ngn_hmedip", "icam_hmedip", "moe_ac3_hmedip", "moe_hmedip")
  foreach(name = names) %dopar% {
    cat(name, "\n")
    cat("Loading hmedips for", name, " \n")
    curr.dset <- loadMedips(dsetToPath(name), 100)    
    cat("Saving dipdata", name,"\n")
    saveMedipsForDipData(curr.dset, name)
  }
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
