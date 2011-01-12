dsetToPath <- function(dsetname) {
  return(paste("data", paste(dsetname, "bed", sep="."), sep="/"))
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

chrToNum <- function(chr) {
  if (chr == "chrX") 
    return(20.0)
  if (chr == "chrY")
    return(21.0)
  if (chr == "chrM")
    return(22.0)
  else
    return(as.double(strsplit(chr, "chr")[[1]][2]))
}

numToChr <- function(chr) {
  if (chr == 20.0) 
    return("chrX")
  if (chr == 21.0)
    return("chrY")
  if (chr == 22.0)
    return("chrM")
  else 
    return(paste("chr", chr, sep=""))
}

chrVecToNum <- function(chr.vec) {
  num.vec <- vector(length=length(chr.vec), mode="double")
  chrs <- c(paste("chr", 1:19, sep=""), "chrX", "chrY", "chrM")
  for (chr in chrs) {
    num.vec[chr.vec == chr] <- chrToNum(chr)
  }
  return(num.vec)
}

# from Software for "Genome-scale DNA methylation mapping of clinical samples at single-nucleotide resolution",
# Gui et al., Nature Methods 2010
fast.fisher = function (x, y = NULL, workspace = 2e+05, hybrid = FALSE, control = list(), 
    or = 1, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95, 
    simulate.p.value = FALSE, B = 2000, cache=F) 
{
    if (nrow(x)!=2 | ncol(x)!=2) stop("Incorrect input format for fast.fisher")
    if (cache) {
      key = paste(x,collapse="_")
      cachedResult = hashTable[[key]]
      if (!is.null(cachedResult)) {
        return(cachedResult)
      }
    }
    # ---- START: cut version of fisher.test ----
    DNAME <- deparse(substitute(x))
    METHOD <- "Fisher's Exact Test for Count Data"
    nr <- nrow(x)
    nc <- ncol(x)
    PVAL <- NULL
    if ((nr == 2) && (nc == 2)) {
        m <- sum(x[, 1])
        n <- sum(x[, 2])
        k <- sum(x[1, ])
        x <- x[1, 1]
        lo <- max(0, k - n)
        hi <- min(k, m)
        NVAL <- or
        names(NVAL) <- "odds ratio"
        support <- lo:hi
        logdc <- dhyper(support, m, n, k, log = TRUE)
        dnhyper <- function(ncp) {
            d <- logdc + log(ncp) * support
            d <- exp(d - max(d))
            d/sum(d)
        }
        mnhyper <- function(ncp) {
            if (ncp == 0) 
                return(lo)
            if (ncp == Inf) 
                return(hi)
            sum(support * dnhyper(ncp))
        }
        pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
            if (ncp == 1) {
                if (upper.tail) 
                  return(phyper(x - 1, m, n, k, lower.tail = FALSE))
                else return(phyper(x, m, n, k))
            }
            if (ncp == 0) {
                if (upper.tail) 
                  return(as.numeric(q <= lo))
                else return(as.numeric(q >= lo))
            }
            if (ncp == Inf) {
                if (upper.tail) 
                  return(as.numeric(q <= hi))
                else return(as.numeric(q >= hi))
            }
            d <- dnhyper(ncp)
            if (upper.tail) 
                sum(d[support >= q])
            else sum(d[support <= q])
        }
        if (is.null(PVAL)) {
            PVAL <- switch(alternative, less = pnhyper(x, or), 
                greater = pnhyper(x, or, upper.tail = TRUE), 
                two.sided = {
                  if (or == 0) 
                    as.numeric(x == lo)
                  else if (or == Inf) 
                    as.numeric(x == hi)
                  else {
                    relErr <- 1 + 10^(-7)
                    d <- dnhyper(or)
                    sum(d[d <= d[x - lo + 1] * relErr])
                  }
                })
            RVAL <- list(p.value = PVAL)
        }
        mle <- function(x) {
            if (x == lo) 
                return(0)
            if (x == hi) 
                return(Inf)
            mu <- mnhyper(1)
            if (mu > x) 
                uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
            else if (mu < x) 
                1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 
                  1))$root
            else 1
        }
        ESTIMATE <- mle(x)
        #names(ESTIMATE) <- "odds ratio"
        RVAL <- c(RVAL, estimate = ESTIMATE, null.value = NVAL)
    }
    RVAL <- c(RVAL, alternative = alternative, method = METHOD, data.name = DNAME)
    attr(RVAL, "class") <- "htest"
    # ---- END: cut version of fisher.test ----    
    if (cache) hashTable[[key]] <<- RVAL # write to global variable
    return(RVAL)                                                                         
}

