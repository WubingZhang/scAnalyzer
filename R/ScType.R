#' Incorporate ScType predictions to Seurat object
#' @docType methods
#' @name ScType
#' @rdname ScType
#'
#' @param query.obj A seurat object for cell type annotation.
#' @param markers A list with Positive (markers$Positive) and Negative markers (markers$Negative) of cell types.
#'
#' @return A seurat object with ScType results attached to meta.data.
#' @export
ScType <- function(query.obj, markers){
  message(Sys.time(), " ScType cell type annotation")
  set.seed(1)
  es.max = sctype_score(GetAssayData(query.obj, slot = "scale.data"), scaled = TRUE,
                        gs = markers$Positive, gs2 = markers$Negative)
  # es.max$assign <- colnames(es.max)[unlist(apply(es.max, 1, which.max))]
  colnames(es.max) <- paste0("ScType.", colnames(es.max))
  rownames(es.max) <- Cells(query.obj)
  query.obj <- AddMetaData(query.obj, metadata = es.max)
  cL_results = do.call("rbind", lapply(as.character(unique(query.obj$seurat_clusters)), function(cl){
    cells <- Cells(query.obj)[as.character(query.obj$seurat_clusters)==cl]
    tmp <- as.matrix(es.max[match(cells, rownames(es.max)), ])
    es.max.cl = colSums(tmp)
    c(seurat_clusters = as.integer(cl), ncells = length(cells), es.max.cl)
  }))
  cL_results <- as.data.frame(cL_results)
  rownames(cL_results) <- cL_results$seurat_clusters
  cL_results <- cL_results[, -1]
  cL_results$assign <- gsub("ScType.", "", colnames(cL_results))[unlist(apply(cL_results[,-1], 1, which.max))+1]
  cL_results$score <- unlist(apply(cL_results[,-1], 1, max))
  cL_results$assign[cL_results$score < cL_results$ncells/4] <- "Unknown"

  idx <- match(as.character(query.obj$seurat_clusters), rownames(cL_results))
  query.obj$ScType.assign <- cL_results[idx, "assign"]
  return(query.obj)
}

#' sctype_score: calculate ScType scores and assign cell types
#' @docType methods
#' @name sctype_score
#' @rdname sctype_score
#'
#' @param scRNAseqData - input scRNA-seq matrix (rownames - genes, column names - cells),
#' @param scale - indicates whether the matrix is scaled (TRUE by default)
#' @param gs - list of gene sets positively expressed in the cell type
#' @param gs2 - list of gene sets that should not be expressed in the cell type (NULL if not applicable)
#'
#' @description GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
#' @author Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#' @return A data frame with ScType scores for each cell type.

sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){

  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T);
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)

  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }

  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]

  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData

  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }

  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]

  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  }))

  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows

  return(as.data.frame(t(es.max)))
}
