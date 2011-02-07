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
  if (lores) {
    diff.enrich <- as.vector(read.table(paste("diff_uniq_genes/lores", de.fname, sep="/"), header=FALSE)$V1)
  } else {
    diff.enrich <- as.vector(read.table(paste("diff_uniq_genes/hires", de.fname, sep="/"), header=FALSE)$V1)
  }
  data.exp1 <- read.table(paste("gene_expr/cufflinks", exp1.fname, sep="/"), header=FALSE, sep="\t")
  colnames(data.exp1) <- c("gene", "rpkm")
  data.exp2 <- read.table(paste("gene_expr/cufflinks", exp2.fname, sep="/"), header=FALSE, sep="\t")
  colnames(data.exp2) <- c('gene', 'rpkm')
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
  plot(log2(expr$exp1+1), log2(expr$exp2+1), pch='.', main=de$diff.enriched.name, col='blue')
}

foldChangeExprBoxPlot.DiffEnrich <- function(de) {
  expr <- getDiffEnrichedExpr.DiffEnrich(de)
  fc <- foldchange(expr$exp1, expr$exp2)
  boxplot(fc, main=de$diff.enriched.name, ylim=c(-60,40))
}

foldChangeExprHist.DiffEnrich <- function(de) {
  expr <- getDiffEnrichedExpr.DiffEnrich(de)
  fc <- foldchange2logratio(foldchange(expr$exp1, expr$exp2))
  hist(fc, main=de$diff.enriched.name)
}

ngnIcamComp <- function(name1, name2, type, lores=FALSE) {
  return(load.DiffEnrich(name1, name2, type, "ngn_expr", "icam_expr", lores=lores))
}

ompNgnComp <- function(name1, name2, type, lores=FALSE) {
  return(load.DiffEnrich(name1, name2, type, "omp_expr", "ngn_expr", lores=lores))
}

drawPlotsHires <- function() {
  par(mfrow=c(2, 3))
  exprPlot.DiffEnrich(ompNgnComp("omp_hmedip_1kb", "ngn_hmedip_1kb", "prom"))
  exprPlot.DiffEnrich(ompNgnComp("omp_hmedip_1kb", "ngn_hmedip_1kb", "prom"))
  exprPlot.DiffEnrich(ompNgnComp("omp_medip_1kb", "ngn_medip_1kb", "prom"))
  
  exprPlot.DiffEnrich(ngnIcamComp("ngn_hmedip_1kb", "icam_hmedip_1kb", "prom"))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_hmedip_1kb", "icam_hmedip_1kb", "gene"))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_medip_1kb", "icam_medip_1kb", "prom"))
}

drawMePlotsLores <- function() {
  par(mfrow=c(2,5))
  exprPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "prom", lores=TRUE))
  exprPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "genes", lores=TRUE))
  exprPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "exons", lores=TRUE))
  exprPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "introns", lores=TRUE))
  exprPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "tss", lores=TRUE))
  
  exprPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "prom", lores=TRUE))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "genes", lores=TRUE))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "exons", lores=TRUE))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "introns", lores=TRUE))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "tss", lores=TRUE))

}

drawHmePlotsLores <- function() {
  par(mfrow=c(2,5))
  exprPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "prom", lores=TRUE))
  exprPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "genes", lores=TRUE))
  exprPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "exons", lores=TRUE))
  exprPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "introns", lores=TRUE))
  exprPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "tss", lores=TRUE))
  
  exprPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "prom", lores=TRUE))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "genes", lores=TRUE))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "exons", lores=TRUE))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "introns", lores=TRUE))
  exprPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "tss", lores=TRUE))
}

drawMeBoxPlotsLores <- function() {
  par(mfrow=c(2,5))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "prom", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "genes", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "exons", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "introns", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_medip", "ngn_medip", "tss", lores=TRUE))
  
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "prom", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "genes", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "exons", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "introns", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_medip", "icam_medip", "tss", lores=TRUE))

}

drawHmeBoxPlotsHires <- function() {
  par(mfrow=c(2,5))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_hmedip_1kb", "ngn_hmedip_1kb", "prom", lores=FALSE))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_hmedip_1kb", "ngn_hmedip_1kb", "gene", lores=FALSE))
  
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_hmedip_1kb", "icam_hmedip_1kb", "prom", lores=FALSE))
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_hmedip_1kb", "icam_hmedip_1kb", "gene", lores=FALSE))
}

drawHmeBoxPlotsLores <- function() {
  par(mfrow=c(2,5))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "prom", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "genes", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "exons", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "introns", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "tss", lores=TRUE))
  
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "prom", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "genes", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "exons", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "introns", lores=TRUE))
  foldChangeExprBoxPlot.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "tss", lores=TRUE))
}

drawHists <- function() {
  par(mfrow=c(4, 2))
  foldChangeExprHist.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "prom", lores=TRUE))
  foldChangeExprHist.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "genes", lores=TRUE))
  foldChangeExprHist.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "exons", lores=TRUE))
  foldChangeExprHist.DiffEnrich(ompNgnComp("omp_hmedip", "ngn_hmedip", "introns", lores=TRUE))
  foldChangeExprHist.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "prom", lores=TRUE))
  foldChangeExprHist.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "genes", lores=TRUE))
  foldChangeExprHist.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "exons", lores=TRUE))
  foldChangeExprHist.DiffEnrich(ngnIcamComp("ngn_hmedip", "icam_hmedip", "introns", lores=TRUE))
}


plotAllExprComparison <- function() {
   icam.expr <- read.table("~/experiment/experiment/stavros_data/gene_expr/cufflinks/icam_expr.txt",sep="\t")
   ngn.expr <- read.table("~/experiment/experiment/stavros_data/gene_expr/cufflinks/ngn_expr.txt", sep="\t")
   omp.expr <- read.table("~/experiment/experiment/stavros_data/gene_expr/cufflinks/omp_expr.txt", sep="\t")
   par(mfrow=c(1,2))
   plot(log2(omp.expr$V2+1), log2(ngn.expr$V2+1), pch='.', col='blue')
   plot(log2(ngn.expr$V2+1), log2(icam.expr$V2+1), pch='.', col='blue')
 }

boxPlotAllExprComparison <- function() {
  icam.expr <- read.table("~/experiment/experiment/stavros_data/gene_expr/cufflinks/icam_expr.txt",sep="\t")
  ngn.expr <- read.table("~/experiment/experiment/stavros_data/gene_expr/cufflinks/ngn_expr.txt", sep="\t")
  omp.expr <- read.table("~/experiment/experiment/stavros_data/gene_expr/cufflinks/omp_expr.txt", sep="\t")
  ngn.icam.fc <- foldchange2logratio(foldchange(ngn.expr$V2, icam.expr$V2))
  omp.ngn.fc <- foldchange2logratio(foldchange(omp.expr$V2, ngn.expr$V2))
  par(mfrow=c(1,2))
  boxplot(omp.ngn.fc, ylim=c(-60,40))
  boxplot(ngn.icam.fc, ylim=c(-60,40))

}



