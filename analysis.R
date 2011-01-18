
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

source("~/src/LEA/medips.R")
source("~/src/LEA/util.R")
source("~/src/LEA/dipdata.R")

euclidDistPlot <- function() {
  pdf("euclid_dist.pdf")
  omp.m <- load.DipData("omp_medip")
  omp.h <- load.DipData('omp_hmedip')
  ngn.m <- load.DipData('ngn_medip')
  ngn.h <- load.DipData('ngn_hmedip')
  moe.ac3.m <- load.DipData('moe_ac3_medip')
  moe.ac3.h <- load.DipData('moe_ac3_hmedip')
  icam.m <- load.DipData('icam_medip')
  icam.h <- load.DipData('icam_hmedip')
  moe.h <- load.DipData('moe_hmedip')
  all.m <- list(omp.m, ngn.m, icam.m, moe.ac3.m, omp.h, ngn.h, icam.h, moe.ac3.h, moe.h)
  all = matrix(nrow=length(all.m), ncol=length(all.m))
  for(i in 1:length(all.m)) {
    for(j in 1:length(all.m)) {
      cat(i,j,'\n')
      all[i,j] <- log2(normDist.DipData(all.m[[i]], all.m[[j]]) + 1)
      cat(all[i,j],'\n')
    }
  }
  rownames(all) <- c('omp m', 'ngn m', 'icam m', 'ac3 m', 'omp h', 'ngn h', 'icam h', 'ac3 h', 'moe h')
  colnames(all) <- c('omp m', 'ngn m', 'icam m', 'ac3 m', 'omp h', 'ngn h', 'icam h', 'ac3 h', 'moe h')
  heatmap.2(all, xlab=rownames(all), ylab=rownames(all))
  dev.off()
}

enrichmentHist <- function() {
  omp.m <- load.DipData("omp_medip")
  omp.h <- load.DipData('omp_hmedip')
  ngn.m <- load.DipData('ngn_medip')
  ngn.h <- load.DipData('ngn_hmedip')
  ac3.m <- load.DipData('moe_ac3_medip')
  icam.m <- load.DipData('icam_medip')
  icam.h <- load.DipData('icam_hmedip')
  ac3.h <- load.DipData('moe_ac3_hmedip')
  pdf("enrichment_hist.pdf")
  par(mfrow=c(3, 4))
  hist(log2(getNormEnrich.DipData(omp.m)), main="omp m")
  hist(log2(getNormEnrich.DipData(ngn.m)), main="ngn m")
  hist(log2(getNormEnrich.DipData(icam.m)), main="icam m")
  hist(log2(getNormEnrich.DipData(ac3.m)), main='ac3 m')
  hist(log2(getNormEnrich.DipData(omp.h)), main="omp h")
  hist(log2(getNormEnrich.DipData(ngn.h)), main="ngn h")
  hist(log2(getNormEnrich.DipData(icam.h)), main="icam h")
  hist(log2(getNormEnrich.DipData(ac3.h)), main='ac3 h')
  dev.off()
}

pairwiseDiffHmeEnrich <- function() {
  cat("omp vs ngn\n")
  dd1 <- load.DipData("moe_hmedip")
  dd2 <- load.DipData("ngn_hmedip")
  p.val <- diffEnrichment.DipData(dd1, dd2)
  p.adj <- mt.rawp2adjp(p.val, proc="BH")$adjp
  write.table(p.adj, file="omp_ngn_pval.txt", quote=FALSE)
  rm(p.val)
  rm(p.adj)
  
  cat("ngn vs icam\n")
  dd1 <- load.DipData("ngn_hmedip")
  dd2 <- load.DipData("icam_hmedip")
  p.val <- diffEnrichment.DipData(dd1, dd2)
  p.adj <- mt.rawp2adjp(p.val, proc="BH")$adjp
  write(p.adj, file="ngn_icam_pval.txt", quote=FALSE)
  rm(p.val)
  rm(p.adj)
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
