library(gtools)
source("common.R")

options(error=recover)

# gene.exp.1 is the gene expression from the start condition
# gene.exp.2 is the gene expression from the end condition
# diff.enriched.genes is a list of the differentially enriched genes
DiffEnrich <- function(diff.enriched.genes, gene.exp.1, gene.exp.2, diff.enriched.name, exp1.name, exp2.name) {
  object <- list(genes=diff.enriched.genes,
                 exp1=gene.exp.1,
                 exp2=gene.exp.2,
                 diff.enriched.name=diff.enriched.name,
                 exp1.name=exp1.name,
                 exp2.name=exp2.name)
  class(object) <- "DiffEnrich"
  return(object)
}

load.DiffEnrich <- function(diff1, diff2, type, exp1.name, exp2.name, lores=FALSE) {
  if (lores) {
    diff.enrich.name <- paste(diff1, diff2, type, "diff_lores_intersect_uniq_genes", sep="_")
  }
  else {
    diff.enrich.name <- paste(diff1, diff2, type, "diff_intersect_uniq_genes", sep="_")
  }
  de.fname <- paste(diff.enrich.name, "txt", sep=".")
  exp1.fname <- paste(exp1.name, "txt", sep=".")
  exp2.fname <- paste(exp2.name, "txt", sep=".")
  old.wd <- getwd()
  setwd(DATA.PATH)
  diff.enrich <- as.vector(read.table(paste("diff_uniq_genes", de.fname, sep="/"), header=FALSE)$V1)
  data.exp1 <- read.table(paste("gene_expr", exp1.fname, sep="/"), header=FALSE, sep="\t")
  colnames(data.exp1) <- c("gene", "rpkm")
  data.exp2 <- read.table(paste("gene_expr", exp2.fname, sep="/"), header=FALSE, sep="\t")
  colnames(data.exp2) <- c("gene", "rpkm")
  obj <- DiffEnrich(diff.enrich, data.exp1, data.exp2, diff.enrich.name, exp1.name, exp2.name)
  setwd(old.wd)
  return(obj)
}

getDiffEnrichedExpr.DiffEnrich <- function(de) {
  genes <- de$genes[which(de$genes %in% de$exp1$gene)]
  exp1 <- de$exp1[which(de$genes %in% de$exp1$gene), "rpkm"]
  exp2 <- de$exp2[which(de$genes %in% de$exp2$gene), "rpkm"]
  out <- data.frame(genes, exp1, exp2)
  colnames(out) <- c("gene", "exp1", "exp2")
  return(out)
}

exprPlot.DiffEnrich <- function(de) {
  expr <- getDiffEnrichedExpr.DiffEnrich(de)
  plot(log2(expr$exp1), log2(expr$exp2), pch=20, main=de$diff.enriched.name)
}

foldChangeExprBoxPlot.DiffEnrich <- function(de) {
  expr <- getDiffEnrichedExpr.DiffEnrich(de)
  fc <- foldchange(expr$exp1, expr$exp2)
  boxplot(fc, ylim=c(-20,20), main=de$diff.enriched.name)
}

foldChangeExprHist.DiffEnrich <- function(de) {
  expr <- getDiffEnrichedExpr.DiffEnrich(de)
  fc <- foldchange2logratio(foldchange(expr$exp1, expr$exp2))
  hist(fc, main=de$diff.enriched.name)
}


#omp.ngn.medip.prom <- load.DiffEnrich("omp_medip_1kb", "ngn_medip_1kb", "prom", "ompgfp", "ngnhigh")
#omp.ngn.medip.gene <- load.DiffEnrich("omp_medip_1kb", "ngn_medip_1kb", "gene", "ompgfp", "ngnhigh")
#omp.ngn.medip.intron <- load.DiffEnrich("omp_medip_1kb", "ngn_medip_1kb", "introns", "ompgfp", "ngnhigh")
#omp.ngn.medip.exon <- load.DiffEnrich("omp_medip_1kb", "ngn_medip_1kb", "exons", "ompgfp", "ngnhigh")

#ngn.icam.hmedip.prom <- load.DiffEnrich("ngn_hmedip_1kb", "icam_hmedip_1kb", "prom", "ompgfp", "ngnhigh")
#ngn.icam.hmedip.gene <- load.DiffEnrich("ngn_hmedip_1kb", "icam_hmedip_1kb", "gene", "ompgfp", "ngnhigh")
#ngn.icam.hmedip.intron <- load.DiffEnrich("ngn_hmedip_1kb", "icam_hmedip_1kb", "introns", "ompgfp", "ngnhigh")
#ngn.icam.hmedip.exon <- load.DiffEnrich("ngn_hmedip_1kb", "icam_hmedip_1kb", "exons", "ompgfp", "ngnhigh")

omp.ngn.medip.prom.lores <-  load.DiffEnrich("omp_hmedip_1kb", "ngn_hmedip_1kb", "prom", "ompgfp", "ngnhigh", lores=FALSE)
omp.ngn.medip.gene.lores <-  load.DiffEnrich("omp_hmedip_1kb", "ngn_hmedip_1kb", "gene", "ompgfp", "ngnhigh", lores=FALSE)
omp.ngn.medip.exon.lores <-  load.DiffEnrich("omp_hmedip_1kb", "ngn_hmedip_1kb", "exons", "ompgfp", "ngnhigh", lores=FALSE)

ngn.icam.hmedip.prom.lores <- load.DiffEnrich("ngn_hmedip_1kb", "icam_hmedip_1kb", "prom", "ompgfp", "ngnhigh", lores=FALSE)
ngn.icam.hmedip.gene.lores <- load.DiffEnrich("ngn_hmedip_1kb", "icam_hmedip_1kb", "gene", "ompgfp", "ngnhigh", lores=FALSE)
ngn.icam.hmedip.exon.lores <- load.DiffEnrich("ngn_hmedip_1kb", "icam_hmedip_1kb", "exons", "ompgfp", "ngnhigh", lores=FALSE)


par(mfrow=c(2,3))
exprPlot.DiffEnrich(omp.ngn.medip.prom.lores)
exprPlot.DiffEnrich(omp.ngn.medip.gene.lores)
exprPlot.DiffEnrich(omp.ngn.medip.exon.lores)

exprPlot.DiffEnrich(ngn.icam.hmedip.prom.lores)
exprPlot.DiffEnrich(ngn.icam.hmedip.gene.lores)
exprPlot.DiffEnrich(ngn.icam.hmedip.exon.lores)


#exprPlot.DiffEnrich(omp.icam.hmedip.prom.lores)
#exprPlot.DiffEnrich(omp.icam.hmedip.gene.lores)

plot(log2(omp.ngn.medip.prom$exp1$rpkm), log2(omp.ngn.medip.prom$exp2$rpkm), pch=20)
