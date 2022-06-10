#' Automatic annotation of cell types using SingleR
#'
#' @docType methods
#' @name scAnnotator.SingleR
#' @rdname scAnnotator.SingleR
#'
#' @param query.obj A seurat object for cell type annotation.
#' @param group.by A character specifying the cluster column in meta.data.
#' @param ref.data A matrix like object specifying the reference data or a Seurat object.
#' @param ref.ann A vector specifying the cell type annotation of the reference data.
#'
#' @author Wubing Zhang
#'
#' @return A data frame with singleR prediction scores.
#'
#' @examples
#'
#' @import Seurat
#' @export
#'
scAnnotator.SingleR <- function(query.obj, group.by = "seurat_clusters", ref.data = NULL, ref.ann = NULL){
  require("SingleR") || stop("Please install SingleR")
  #### Retrieve reference data ####
  if(!is.null(ref.data)){
    if(class(ref.data)[1] == "Seurat"){
      ref.data <- subset(ref.data, downsample = 200)
      ref.ann <- Idents(ref.data)
      ref.data <- as.matrix(GetAssayData(object = ref.data, slot = "data"))
    }
    ref.data <- as.matrix(ref.data)
  }
  if(is.null(ref.data)){
    require("celldex") || stop("Please install celldex")
    ## Combine multiple references
    hpca <- suppressMessages(celldex::HumanPrimaryCellAtlasData())
    Encode <- suppressMessages(celldex::BlueprintEncodeData())
    Monaco <- suppressMessages(celldex::MonacoImmuneData())
    # hpca$label.main <- paste0("HPCA.", hpca$label.main)
    # Encode$label.main <- paste0("Encode.", Encode$label.main)
    # Monaco$label.main <- paste0("Monaco.", Monaco$label.main)
    shared <- intersect(intersect(rownames(hpca), rownames(Encode)), rownames(Monaco))
    ref.data <- cbind(hpca[shared,], Encode[shared,])
    ref.data <- cbind(ref.data, Monaco[shared,])
    ref.ann <- ref.data$label.main
    ref.ann[grepl("Plasma_cell|Plasma cells", ref.data$label.fine)] <- "Plasma"
    ref.ann[grepl("CD8", ref.ann)] <- "CD8T"
    ref.ann[grepl("CD4", ref.ann)] <- "CD4T"
    ref.ann[grepl("NK", ref.ann)] <- "NK"
    ref.ann[grepl("B cell|B-cell|B_cell", ref.ann)] <- "B"
    ref.ann[grepl("Epithelial", ref.ann)] <- "Epithelial"
    ref.ann[grepl("Endothelial", ref.ann)] <- "Endothelial"
    ref.ann[grepl("Fibroblast", ref.ann)] <- "Fibroblast"
    ref.ann[grepl("Dendritic", ref.ann)] <- "DC"
    ref.ann[grepl("Macrophage", ref.ann)] <- "Macrophage"
    ref.ann[ref.ann %in% names(table(ref.ann))[table(ref.ann)<5]] <- "Other"
    ref.ann[grepl("HSC|muscle|Neuron|stem|Gametocytes|MSC|Chondrocytes|Monocyte|iPS|Platelets|BM|Eryth|Osteoblasts|Adipocytes|T_cells|T cells", ref.ann)] <- "Other"
    idx <- ref.ann!="Other"
    ref.data <- ref.data[, idx]
    ref.ann <- ref.ann[idx]
    ref.data$label.main <- ref.ann
  }
  ## Annotation
  vargenes <- VariableFeatures(query.obj)
  pred.singleR <- suppressWarnings(
    SingleR::SingleR(test = as.matrix(GetAssayData(query.obj, slot = "data"))[vargenes,],
                     clusters = query.obj[[]][[group.by]],
                     ref = ref.data, labels = ref.ann, de.method = "wilcox")
  )
  scores <- pred.singleR[match(query.obj[[]][[group.by]], rownames(pred.singleR)), ]
  rownames(scores) <- colnames(query.obj)
  return(scores)
}
