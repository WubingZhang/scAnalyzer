#' Annotating cell types using SingleR
#'
#' @docType methods
#' @name CellTypeAnnotator
#' @rdname CellTypeAnnotator
#'
#' @param query.obj A seurat object for cell type annotation.
#' @param clusters A character specifying the cluster column in meta.data.
#' @param method A vector specifying the cell type annotation methods.
#' Should be one of singler, sclearn, scmap-cell, scmap-cluster, seurat
#' @param ref.data A matrix like object specifying the reference data or a Seurat object.
#' @param ref.ann A vector specifying the cell type annotation of the reference data.
#' @param plot A boolean indicating whether to save and visualize the prediction results.
#' @param outdir Path to the output directory.
#'
#' @author Wubing Zhang
#'
#' @return A seurat object with prediction results added into meta.data.
#'
#' @examples
#'
#' @import Seurat SingleR celldex scLearn scmap SingleCellExperiment
#' @export
#'
CellTypeAnnotator <- function(query.obj,
                              clusters = "seurat_clusters",
                              method = c("scmap-cell", "seurat"),
                              ref.data = NULL,
                              ref.ann = NULL,
                              plot = TRUE,
                              outdir = "./"){
  if(length(names(ref.ann))==0){ names(ref.ann) <- colnames(ref.data) }
  requireNamespace("Seurat") || stop("Please install Seurat")

  #### Retrieve reference data ####
  ref.obj <- NULL
  if(class(ref.data) == "Seurat"){
    ref.obj <- ref.data
    ref.data <- as.matrix(GetAssayData(object = ref.obj, slot = "data"))
  }
  ref.data <- as.matrix(ref.data)

  #### SingleR Cell type annotation ####
  if("singler" %in% tolower(method)){
    requireNamespace("SingleR") || stop("Please install SingleR")
    message(Sys.time(), " SingleR cell type annotation")
    query.sce <- Seurat::as.SingleCellExperiment(query.obj)
    if(is.null(ref.data)){
      requireNamespace("celldex") || stop("Please install celldex")
      ## Combine multiple references
      hpca <- celldex::HumanPrimaryCellAtlasData()
      Encode <- celldex::BlueprintEncodeData()
      Monaco <- celldex::MonacoImmuneData()
      hpca$label.main <- paste0("HPCA.", hpca$label.main)
      Encode$label.main <- paste0("Encode.", Encode$label.main)
      Monaco$label.main <- paste0("Monaco.", Monaco$label.main)

      shared <- intersect(intersect(rownames(hpca), rownames(Encode)), rownames(Monaco))
      combined <- cbind(hpca[shared,], Encode[shared,])
      combined <- cbind(combined, Monaco[shared,])

      ## Annotation
      pred.singleR <- suppressWarnings(
        SingleR::SingleR(test = query.sce, clusters = query.obj[[]][[clusters]], ref = combined,
                labels = combined$label.main, de.method = "wilcox")
        )
      query.obj[["singleR.assign"]] <- pred.singleR$labels[match(query.obj[[]][[clusters]], rownames(pred.singleR))]
      query.obj[["singleR.score1"]] <- pred.singleR$tuning.scores[match(query.obj[[]][[clusters]], rownames(pred.singleR)), 1]
      query.obj[["singleR.score2"]] <- pred.singleR$tuning.scores[match(query.obj[[]][[clusters]], rownames(pred.singleR)), 2]
      query.obj[["singleR.assign"]][query.obj[["singleR.score1"]]<0.3] = "Other"
    }else{
      pred.singleR <- suppressWarnings(
        SingleR::SingleR(test = query.sce, clusters = query.obj[[]][[clusters]],
                         ref = ref.data, labels = ref.ann, de.method = "wilcox"))
      query.obj[["singleR.assign"]] <- pred.singleR$labels[match(query.obj[[]][[clusters]], rownames(pred.singleR))]
      query.obj[["singleR.score1"]] <- pred.singleR$tuning.scores[match(query.obj[[]][[clusters]], rownames(pred.singleR)), 1]
      query.obj[["singleR.score2"]] <- pred.singleR$tuning.scores[match(query.obj[[]][[clusters]], rownames(pred.singleR)), 2]
      query.obj[["singleR.assign"]][query.obj[["singleR.score1"]]<0.3] = "Other"
    }
    p <- DimPlot(query.obj, reduction = "umap", group.by = "singleR.assign", label=TRUE, label.size = 6)
    if(plot){
      saveRDS(pred.singleR, paste0(outdir, "/SingleR_assignment.rds"))
      ggsave(plot = p, paste0(outdir, "/UmapPlot_SingleR_Assign_", Sys.Date(), ".pdf"), width = 8, height = 6)
    }
  }
  if("sclearn" %in% tolower(method)){
    requireNamespace("scLearn") || stop("Please install scLearn")
    message(Sys.time(), " sclearn cell type annotation")
    # cell quality control and rare cell type filtered and feature selection
    ref_filtered <- scLearn::Cell_type_filter(ref.data, ref.ann, min_cell_number = 10)
    ref_filtered$expression_profile[ref_filtered$expression_profile<0] <- 0
    varGenes <- scLearn::Feature_selection_M3Drop(ref_filtered$expression_profile, log_normalized = TRUE)
    scLearn_mod <- scLearn::scLearn_model_learning(varGenes, ref_filtered$expression_profile,
                                          ref_filtered$sample_information_cellType,
                                          bootstrap_times=10)
    querydata <- GetAssayData(object = query.obj, slot = "data")
    querydata[querydata<0] <- 0
    scLearn_pred <- scLearn::scLearn_cell_assignment(scLearn_mod, querydata, diff=0.05,
                                            threshold_use=TRUE, vote_rate=0.6)
    rownames(scLearn_pred) <- scLearn_pred$Query_cell_id
    query.obj[["scLearn.assign"]] <- scLearn_pred[colnames(query.obj), "Predict_cell_type"]

    p <- DimPlot(query.obj, reduction = "umap", group.by = "scLearn.assign", label=TRUE, label.size = 6)
    if(plot){
      saveRDS(scLearn_pred, paste0(outdir, "/scLearn_assignment.rds"))
      ggsave(plot = p, paste0(outdir, "/UmapPlot_scLearn_Assign_", Sys.Date(), ".pdf"), width = 8, height = 6)
    }
  }
  if("seurat" %in% tolower(method)){
    require(Seurat)
    message(Sys.time(), " seurat cell type annotation")
    if(is.null(ref.obj)){
      ref.obj <- CreateSeuratObject(counts = ref.data, project = "ref")
      ref.obj <- SetAssayData(object = ref.obj, slot = "data", new.data = ref.data)
    }
    ref.obj <- FindVariableFeatures(ref.obj, selection.method = "vst", nfeatures = 2000)
    ref.obj <- ScaleData(ref.obj, verbose = FALSE)
    ref.obj <- RunPCA(ref.obj, npcs = 30, verbose = FALSE)
    anchors <- FindTransferAnchors(reference = ref.obj, query = query.obj,
                                   dims = 1:30, reference.reduction = "pca")
    predictions <- TransferData(anchorset = anchors, refdata = ref.ann, dims = 1:30)
    tmp <- predictions[, c("predicted.id", "prediction.score.max")]
    colnames(tmp) <- c("Seurat.assign", "Seurat.score")
    query.obj <- AddMetaData(query.obj, metadata = tmp)
    # save the prediction results
    p <- DimPlot(query.obj, reduction = "umap", group.by = "Seurat.assign", label=TRUE, label.size = 6)
    if(plot){
      saveRDS(predictions, paste0(outdir, "/Seurat_assignment.rds"))
      ggsave(plot = p, paste0(outdir, "/UmapPlot_Seurat_Assign_", Sys.Date(), ".pdf"), width = 8, height = 6)
    }
  }
  if("scmap-cell" %in% tolower(method)){
    requireNamespace("scmap") || stop("Please install scmap")
    requireNamespace("SingleCellExperiment") || stop("Please install SingleCellExperiment")
    message(Sys.time(), " scmap-cell cell type annotation")
    set.seed(1)
    ref.sce <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = as.matrix(ref.data),
                                                  logcounts = as.matrix(ref.data)),
                                    colData = DataFrame(label=ref.ann),
                                    rowData = DataFrame(feature_symbol=rownames(ref.data)))
    query.sce <- Seurat::as.SingleCellExperiment(query.obj)
    rowData(query.sce) <- DataFrame(feature_symbol=rownames(query.obj))

    # feature_symbol column in the rowData slot
    ref.sce <- suppressMessages(suppressWarnings(scmap::selectFeatures(ref.sce, suppress_plot = TRUE)))
    ref.sce <- suppressMessages(suppressWarnings(scmap::indexCell(ref.sce)))
    scmapCell_results <- scmap::scmapCell(query.sce, list(ref = metadata(ref.sce)$scmap_cell_index))
    scmapCell_clusters <- scmap::scmapCell2Cluster(scmapCell_results,
      list(as.character(colData(ref.sce)$label))
    )
    tmp <- data.frame(scmap.cell.assign = scmapCell_clusters$scmap_cluster_labs,
                      scmap.cell.score = scmapCell_clusters$scmap_cluster_siml,
                      row.names = colnames(query.sce))
    colnames(tmp) <- c("scmapCell.assign", "scmapCell.score")
    query.obj <- AddMetaData(query.obj, metadata = tmp)
    # save the prediction results
    p <- DimPlot(query.obj, reduction = "umap", group.by = "scmapCell.assign", label=TRUE, label.size = 6)
    if(plot){
      saveRDS(scmapCell_clusters, paste0(outdir, "/scmapCell_assignment.rds"))
      ggsave(plot = p, paste0(outdir, "/UmapPlot_scmapCell_Assign_", Sys.Date(), ".pdf"), width = 8, height = 6)
    }
  }
  if("scmap-cluster" %in% tolower(method)){
    requireNamespace("scmap") || stop("Please install scmap")
    requireNamespace("SingleCellExperiment") || stop("Please install SingleCellExperiment")
    message(Sys.time(), " scmap-cluster cell type annotation")
    set.seed(1)
    ref.sce <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = as.matrix(ref.data),
                                                  logcounts = as.matrix(ref.data)),
                                    colData = DataFrame(label=ref.ann),
                                    rowData = DataFrame(feature_symbol=rownames(ref.data)))
    query.sce <- as.SingleCellExperiment(query.obj)
    rowData(query.sce) <- DataFrame(feature_symbol=rownames(query.obj))

    ref.sce <- suppressMessages(suppressWarnings(scmap::selectFeatures(ref.sce, suppress_plot = TRUE)))
    ref.sce <- suppressMessages(suppressWarnings(scmap::indexCluster(ref.sce, cluster_col = "label")))
    scmapCluster_results <- scmap::scmapCluster(query.sce, list(ref = metadata(ref.sce)$scmap_cluster_index))

    tmp <- data.frame(scmap.cluster.assign = scmapCluster_results$scmap_cluster_labs,
                      scmap.cluster.score = scmapCluster_results$scmap_cluster_siml,
                      row.names = colnames(query.sce))
    colnames(tmp) <- c("scmapCluster.assign", "scmapCluster.score")
    query.obj <- AddMetaData(query.obj, metadata = tmp)
    # save the prediction results
    p <- DimPlot(query.obj, reduction = "umap", group.by = "scmapCluster.assign", label=TRUE, label.size = 6)
    if(plot){
      saveRDS(scmapCluster_results, paste0(outdir, "/scmapCluster_assignment.rds"))
      ggsave(plot = p, paste0(outdir, "/UmapPlot_scmapCluster_Assign_", Sys.Date(), ".pdf"), width = 8, height = 6)
    }
  }
  return(query.obj)
}
