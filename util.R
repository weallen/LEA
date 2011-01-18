callSignificant <- function(p.vals, fdr.level=0.1) {
  qvals <- try(qvalue(signif(p.vals, 8)))
  if ("qvalues" %in% names(qvals)) {
    cat("using qvals\n")
    return(qvals$qvalues)
  } else {
    cat("Doing Bonferroni correction\n")
    return(pmin(p.vals * length(p.vals), 1))
  }
}

# from Storey's qvalue package, with gui stuff removed
qvalue <- function(p=NULL, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=NULL, robust=FALSE, smooth.df = 3, smooth.log.pi0 = FALSE) {
#Input
#=============================================================================
#p: a vector of p-values (only necessary input)
#fdr.level: a level at which to control the FDR (optional)
#lambda: the value of the tuning parameter to estimate pi0 (optional)
#pi0.method: either "smoother" or "bootstrap"; the method for automatically
#           choosing tuning parameter in the estimation of pi0, the proportion
#           of true null hypotheses
#robust: an indicator of whether it is desired to make the estimate more robust
#        for small p-values and a direct finite sample estimate of pFDR (optional)
#smooth.df: degrees of freedom to use in smoother (optional)
#smooth.log.pi0: should smoothing be done on log scale? (optional)
#
#Output
#=============================================================================
#call: gives the function call
#pi0: an estimate of the proportion of null p-values
#qvalues: a vector of the estimated q-values (the main quantity of interest)
#pvalues: a vector of the original p-values
#significant: if fdr.level is specified, an indicator of whether the q-value
#    fell below fdr.level (taking all such q-values to be significant controls
#    FDR at level fdr.level)


#This is just some pre-processing
    if(min(p)<0 || max(p)>1) {
      print("ERROR: p-values not in valid range.")
      return(0)
    }
    if(length(lambda)>1 && length(lambda)<4) {      
      print("ERROR: If length of lambda greater than 1, you need at least 4 values.")
      return(0)
    }
    if(length(lambda)>1 && (min(lambda) < 0 || max(lambda) >= 1)) { ## change by Alan:  check for valid range for lambda
      print("ERROR: Lambda must be within [0, 1).")
      return(0)
    }
    m <- length(p)
#These next few functions are the various ways to estimate pi0
    if(length(lambda)==1) {
        if(lambda<0 || lambda>=1) { ## change by Alan:  check for valid range for lambda          
          print("ERROR: Lambda must be within [0, 1).")
          return(0)
        }

        pi0 <- mean(p >= lambda)/(1-lambda)
        pi0 <- min(pi0,1)
    }
    else {
        pi0 <- rep(0,length(lambda))
        for(i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])
        }

        if(pi0.method=="smoother") {
            if(smooth.log.pi0)
              pi0 <- log(pi0)

            spi0 <- smooth.spline(lambda,pi0,df=smooth.df)
            pi0 <- predict(spi0,x=max(lambda))$y

            if(smooth.log.pi0)
              pi0 <- exp(pi0)
            pi0 <- min(pi0,1)
        }
        else if(pi0.method=="bootstrap") {
            minpi0 <- min(pi0)
            mse <- rep(0,length(lambda))
            pi0.boot <- rep(0,length(lambda))
            for(i in 1:100) {
                p.boot <- sample(p,size=m,replace=TRUE)
                for(i in 1:length(lambda)) {
                    pi0.boot[i] <- mean(p.boot>lambda[i])/(1-lambda[i])
                }
                mse <- mse + (pi0.boot-minpi0)^2
            }
            pi0 <- min(pi0[mse==min(mse)])
            pi0 <- min(pi0,1)
        }
        else {  ## change by Alan: check for valid choice of 'pi0.method' (only necessary on command line)
            print("ERROR: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
            return(0)
        }
    }
    if(pi0 <= 0) {
      print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
      return(0)
    }
    if(!is.null(fdr.level) && (fdr.level<=0 || fdr.level>1)) {  ## change by Alan:  check for valid fdr.level      
      print("ERROR: 'fdr.level' must be within (0, 1].")
      return(0)
    }
#The estimated q-values calculated here
    u <- order(p)

    # change by Alan
    # ranking function which returns number of observations less than or equal
    qvalue.rank <- function(x) {
      idx <- sort.list(x)

      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
 
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl

      return(tbl)
    }

    v <- qvalue.rank(p)
    
    qvalue <- pi0*m*p/v
    if(robust) {
        qvalue <- pi0*m*p/(v*(1-(1-p)^m))
    }
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
    }
#The results are returned
    if(!is.null(fdr.level)) {
        retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, fdr.level=fdr.level, ## change by Alan
          significant=(qvalue <= fdr.level), lambda=lambda)
    }
    else {
        retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, lambda=lambda)
    }
    class(retval) <- "qvalue"
    return(retval)
}


dsetToPath <- function(dsetname) {
  return(paste("/gpfs/home/wallen/data/aligned_reads/", paste(dsetname, "bed", sep="."), sep="/"))
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

