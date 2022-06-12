#' Normalization of single-cell RNA-seq read counts
#'
#' @docType methods
#' @name scNormalize
#' @rdname scNormalize
#'
#' @param counts A matrix.
#' @param norm.method One of LogNormalize, SCTransform, downsample, quminorm, and census.
#' @param mc.cores An integer, specifying the number of cores for quminorm and downsample normalization.
#' @param fraction A numeric between 0 and 1, specifying the fraction of cells to be downsampling normalized.
#'
#' @author Wubing Zhang
#'
#' @return A normalized count matrix.
#'
#' @examples
#'
#' @export
#'
scNormalize <- function(counts,
                        norm.method = c("LogNormalize", "SCT",
                                        "downsample",
                                        "quminorm",
                                        "census"),
                        mc.cores=4,
                        fraction = 0.85){
  counts[is.na(counts)] <- 0
  if("SCT" %in% norm.method){
    normcounts <- sctransform::vst(counts, return_corrected_umi = TRUE,
                                   verbosity = FALSE)$umi_corrected
  }
  if("LogNormalize" %in% norm.method){
    normcounts <- Seurat::LogNormalize(counts, verbose = FALSE)
    normcounts <- exp(normcounts)-1
  }
  if("quminorm" %in% norm.method){
    normcounts <- quminorm::quminorm(counts, mc.cores=mc.cores)
  }
  if("downsample" %in% norm.method){
    set.seed(42)
    libSizes <- colSums(counts)
    targetSize <- floor(quantile(libSizes, 1-fraction))
    toRemove <- libSizes-targetSize
    ii <- which(toRemove>0)
    normcounts <- counts
    for(i in ii){
      readsGet <- sort(sample(seq_len(sum(counts[,i])), toRemove[i]))
      cumCounts <-  c(0, cumsum(counts[, i]))
      tmp <- hist(readsGet, breaks = cumCounts, plot=FALSE)$count
      normcounts[,i] <- normcounts[,i] - tmp
    }
  }
  if("census" %in% norm.method){
    normcounts <- t(t(counts) / colSums(counts>0)) * median(colSums(counts>0))
  }
  return(normcounts)
}
