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
      gene_order_file <- system.file("extdata", "geneorder_hg20.txt.gz", package="scAnalyzer")
      gene_order_file <- gzfile(gene_order_file, "r")
    }
    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = as.matrix(SeuratObj@assays$RNA@counts),
                                         gene_order_file = gene_order_file,
                                         annotations_file=paste0(outdir, "/annotation_for_infercnv.txt"),
                                         ref_group_names = normal_groups)
    infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=outdir, cluster_by_groups=TRUE, denoise=TRUE, HMM=TRUE)
    SeuratObj <- add_to_seurat(SeuratObj, outdir, top_n = 10, bp_tolerance = 2000000)

    ## Visualize the infercnv prediction

  }
  if("numbat" %in% tolower(methods)){
    ## Add here
  }
  return(SeuratObj)
}
