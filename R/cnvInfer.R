#' Inferring CNV using different methods
#'
#' @docType methods
#' @name cnvInfer
#' @rdname cnvInfer
#'
#' @param SeuratObj A seurat object including count matrix and cell cluster annotation.
#' @param ann.column An integer or a character specifying the cell cluster column in the seurat meta.data.
#' @param normal_groups A vector specifying the groups of normal cells, which mush match the name shown in meta.data.
#' @param methods copykat, infercnv, numbat
#' @param outdir path to output intermediate results.
#' @param gene_order_file path to the gene order file (required for infercnv); using geneorder_hg20 by default.
#'
#' @author Wubing Zhang
#'
#' @return A seurat object with prediction results added into meta.data.
#'
#' @examples
#'
#' @import Seurat infercnv copykat
#' @export
cnvInfer <- function(SeuratObj,
                     ann.column = "seurat_clusters",
                     normal_groups = NULL,
                     methods = c("copykat", "infercnv"),
                     outdir = "./",
                     gene_order_file = NULL){

  if("copykat" %in% tolower(methods)){
    require(copykat)
    if(is.null(normal_groups) | is.null(ann.column)){
      norm.cells <- subset(x = SeuratObj, subset = PTPRC > 3) %>% colnames()
    }else{
      norm.cells <- rownames(SeuratObj@meta.data)[SeuratObj@meta.data[,ann.column] %in% normal_groups]
    }
    if(length(norm.cells)<10){
      warning("Too few specified normal cells !")
      copykat_res <- copykat(rawmat=as.matrix(SeuratObj@assays$RNA@counts),
                             win.size=25, sam.name=outdir, output.seg="FLASE",
                             plot.genes="TRUE", genome="hg20", n.cores = 16)
    }else{
      copykat_res <- copykat(rawmat=as.matrix(SeuratObj@assays$RNA@counts),
                             win.size=25, sam.name=outdir,
                             norm.cell.names=norm.cells, output.seg="FLASE",
                             plot.genes="TRUE", genome="hg20", n.cores = 16)
    }
    pred <- data.frame(copykat_res$prediction)
    SeuratObj@meta.data$copykat = pred[rownames(SeuratObj@meta.data), "copykat.pred"]

    ## Visualize the copykat prediction
  }

  if("infercnv" %in% tolower(methods)){
    require(infercnv)
    ann = SeuratObj@meta.data
    ann = cbind(rownames(ann), ann[,ann.column])
    write.table(ann, paste0(outdir, "/annotation_for_infercnv.txt"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    if(is.null(gene_order_file)){
      gene_order_file <- system.file("extdata", "geneorder_hg20.txt", package="scAnalyzer")
      gene_order_file <- gzipfile(gene_order_file, "r")
    }
    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = as.matrix(SeuratObj@assays$RNA@counts),
                                         gene_order_file = gene_order_file,
                                         annotations_file=paste0(outdir, "/annotation_for_infercnv.txt"),
                                         ref_group_names = normal_groups)
    infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=outdir, cluster_by_groups=TRUE, denoise=TRUE, HMM=TRUE)

    require(dplyr)
    require(EnvStats)

    ref_file <- paste0(out_dir, '/infercnv.references.txt')
    obs_file <- paste0(out_dir, '/infercnv.observations.txt')
    ref <- read.csv(ref_file, header=T, sep = ' ', row.names = 1)
    obs <- read.csv(obs_file, header=T, sep = ' ', row.names = 1)

    ref_var <- mean(as.numeric(ref %>% summarise_if(is.numeric, var)))

    df <- data.frame()
    for (cell in colnames(observation)) {
      observation_var <- var(observation[, cell])
      log_p <- pf(observation_var / reference_var, nrow(observation) - 1, nrow(reference) - 1, lower.tail = FALSE, log.p = TRUE)

      # Bonferroni adjusted p-value
      adjusted_log_p <- min(log_p + log(ncol(observation)), 0)
      if (adjusted_log_p <= log(0.05)) {
        type <- "T"
      } else {
        type <- "N"
      }

      df <- rbind(df, data.frame("cell" = cell, "var" = observation_var, "negative_log_p" = -adjusted_log_p, "type" = type))
    }
    row.names(df) <- df[["cell"]]
    SeuratObj@meta.data$infercnv = df[rownames(SeuratObj@meta.data), "type"]
    SeuratObj@meta.data$infercnv[is.na(SeuratObj@meta.data$infercnv)] <- "N"

    ## Visualize the infercnv prediction
    DimPlot(SeuratObj, group.by = "infercnv")

  }
  if("numbat" %in% tolower(methods)){
    ## Add here
  }
  return(SeuratObj)
}
