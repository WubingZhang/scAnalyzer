#' Inferring CNV using different methods
#'
#' @docType methods
#' @name infercnvParser
#' @rdname infercnvParser
#'
#' @param infercnv_obj An Infercnv object.
#' @param window.width width of windows to group genes
#'
#' @author Wubing Zhang
#'
#' @return A data frame with inferred malignancy.
#'
#' @examples
#'
#' @export
infercnvParser <- function(infercnv_obj, window.width = 1){
  expr.data <- cbind(infercnv_obj@gene_order, infercnv_obj@expr.data)
  if(window.width>100){
    expr.data$start <- floor(expr.data$start / window.width)
    expr.data <- aggregate(expr.data[,-(1:3)], by = list(chr = expr.data$chr, start = expr.data$start), mean)
    rownames(expr.data) <- paste0(expr.data$chr, "_", expr.data$start)
  }
  expr.data <- expr.data[,-(1:2)]
  ref_var <- mean(apply(expr.data[, unlist(infercnv_obj@reference_grouped_cell_indices)], 2, var))
  pvals <- apply(expr.data, 2, function(x){
    tmp <- varTest(x, alternative = "greater", sigma.squared = ref_var)
    pchisq(tmp$statistic, tmp$parameters, ncp = 0, lower.tail = FALSE, log.p = FALSE)
  })
  Padj <- p.adjust(pvals, method = "BH")
  results <- data.frame(CID = colnames(expr.data), InferCNV.Pval = pvals, InferCNV.FDR = Padj)
  results$Infer.Tumor = "N"
  results$Infer.Tumor[results$InferCNV.FDR<0.01] = "T"
  return(results)
}
