library(biomaRt)

loadMouseGenes <- function() {
  if (!file.exists("genes.biomart")) {
    ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    mm9.genes <- getBM(attributes = c("chromosome_name", "strand", "start_position",
                         "end_position", "description"), mart=ensembl)
    szve(mm9.genes, file="genes.biomart")
  } else {
    load("genes.biomart")
  }
  return(mm9.genes)  
}

overlap <- function(pos, trans.start, dist=1000) {
  if ((pos > trans.start - dist) && (pos < trans.start + dist)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

