#' Automatic annotation of cell types using Seurat label transfer
#'
#' @docType methods
#' @name scAnnotator.Seurat
#' @rdname scAnnotator.Seurat
#'
#' @param query.obj A seurat object for querying
#' @param ref.obj A seurat object as reference
#' @param norm.method SCT or LogNormalize
#'
#' @author Wubing Zhang
#'
#' @return A data frame with Seurat prediction scores.
#'
#' @examples
#'
#' @export
#'
scAnnotator.Seurat <- function(query.obj, ref.obj, norm.method = "SCT"){
  require(Seurat)
  if(norm.method=="SCT"){
    if(DefaultAssay(ref.obj)!="SCT")
      ref.obj <- ref.obj %>% SCTransform(verbose = FALSE)
    if(DefaultAssay(query.obj)!="SCT")
      query.obj <- query.obj %>% SCTransform(verbose = FALSE)
  }else{
    ref.obj <- ref.obj %>% NormalizeData(verbose = FALSE) %>%
      FindVariableFeatures(verbose = FALSE) %>% ScaleData(verbose = FALSE)
    query.obj <- query.obj %>% NormalizeData(verbose = FALSE) %>%
      FindVariableFeatures(verbose = FALSE) %>% ScaleData(verbose = FALSE)
  }
  ref.obj <- RunPCA(ref.obj, npcs = 30, verbose = FALSE)
  anchors <- FindTransferAnchors(reference = ref.obj, query = query.obj,
                                 normalization.method = "SCT",
                                 reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors,
                              refdata = Idents(ref.obj), dims = 1:30)
  return(predictions)
}
