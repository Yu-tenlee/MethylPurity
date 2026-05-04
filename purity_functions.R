#' Calculate Tumor Purity based on DNA Methylation Data
#'
#' This package provides functions to estimate tumor purity using
#' differentially methylated CpG sites (iDMCs).
#'
#' @keywords internal
"_PACKAGE"

#' Calculate Variance
#' @param x A numeric vector.
#' @return Variance of x.
myVar <- function(x){
  var(x, na.rm = TRUE)
}

#' Wilcoxon Rank Sum Test for iDMC selection
#' @param x A numeric vector of methylation values.
#' @param tumor_idx Column indices for tumor samples.
#' @param normal_idx Column indices for normal samples.
#' @return P-value of the test.
myRanksum <- function(x, tumor_idx, normal_idx){
  if (length(na.omit(x[tumor_idx])) == 0 | length(na.omit(x[normal_idx])) == 0){
    return(NA)
  } else {
    pval <- suppressWarnings(wilcox.test(as.numeric(x[tumor_idx]), as.numeric(x[normal_idx]))$p.value)
    return(pval)
  }
}

#' Get the peak value from density distribution
#' @param dat A numeric vector.
#' @return The x value corresponding to the highest density peak.
get_peak <- function(dat){
  d <- density(dat, na.rm = TRUE, kernel = "gaussian")
  d$x[which.max(d$y)]
}

#' Identify intrinsic Differentially Methylated CpGs (iDMC)
#' @param tumor.data Matrix of tumor methylation data.
#' @param normal.data Matrix of normal methylation data.
#' @return A character vector of iDMC probe names.
get_iDMC <- function(tumor.data, normal.data){
  all.dat <- cbind(tumor.data, normal.data)
  tumor.sample <- colnames(tumor.data)
  normal.sample <- colnames(normal.data)

  # Use indices for faster processing
  t_idx <- 1:ncol(tumor.data)
  n_idx <- (ncol(tumor.data) + 1):ncol(all.dat)

  ranksum.pval <- apply(all.dat, 1, function(x) myRanksum(x, t_idx, n_idx))
  tumor.var <- apply(tumor.data, 1, myVar)

  out <- data.frame(ranksum.pval, tumor.var)
  cDMC <- out[!is.na(out$ranksum.pval) & out$tumor.var >= 0.005, ]

  if(nrow(cDMC) == 0) return(character(0)) # Handle empty case

  idx <- order(cDMC$ranksum.pval, decreasing = FALSE)[1:min(1000, nrow(cDMC))]
  iDMC <- rownames(cDMC)[idx]
  return(iDMC)
}

#' Main Function: Calculate Tumor Purity
#' @param tumor.data A matrix of tumor beta values (rows: probes, cols: samples).
#' @param normal.data A matrix of normal beta values (rows: probes, cols: samples).
#' @return A named numeric vector of purity values.
#' @export
getPurity <- function(tumor.data, normal.data){
  iDMC <- get_iDMC(tumor.data, normal.data)

  if(length(iDMC) == 0) {
    stop("No iDMCs found. Try adjusting variance threshold.")
  }

  tumor.sample <- colnames(tumor.data)
  normal.sample <- colnames(normal.data)

  idmc.dat <- data.frame(tumor.data[iDMC, , drop=FALSE], normal.data[iDMC, , drop=FALSE])
  colnames(idmc.dat) <- c(tumor.sample, normal.sample)

  idmc.dat$hyper <- rowMeans(idmc.dat[, tumor.sample, drop=FALSE], na.rm = TRUE) >
    rowMeans(idmc.dat[, normal.sample, drop=FALSE], na.rm = TRUE)

  cat("Calculating tumor purity ...\n")
  purity <- rep(NA, length(tumor.sample))
  names(purity) <- tumor.sample

  for(t in tumor.sample){
    beta.adj <- c(idmc.dat[idmc.dat$hyper == TRUE, t],
                  1 - idmc.dat[idmc.dat$hyper == FALSE, t])
    pu <- get_peak(beta.adj)
    purity[t] <- pu
  }
  return(purity)
}
