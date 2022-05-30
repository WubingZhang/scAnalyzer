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
#' @param nthreads An integer specifying the number of threads for running inferCNV
#' @param outdir path to output intermediate results.
#' @param gene_order_file path to the gene order file (required for infercnv); using geneorder_hg20 by default.
#'
#' @author Wubing Zhang
#'
#' @return A seurat object with prediction results added into meta.data.
#'
#' @examples
#'
#' @import Seurat infercnv
#' @export
cnvInfer <- function(SeuratObj,
                     ann.column = "seurat_clusters",
                     normal_groups = NULL,
                     methods = c("infercnv"),
                     nthreads = 8,
                     outdir = "./",
                     gene_order_file = NULL){
  if("copykat" %in% tolower(methods)){
    require(copykat)
    if(is.null(normal_groups) | is.null(ann.column)){
      norm.cells <- subset(x = SeuratObj, subset = EPCAM < 0.5) %>% colnames()
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
    require(dplyr)
    require(EnvStats)

    SeuratObj$Clusters <- as.character(SeuratObj@meta.data[, ann.column])
    ann = data.frame(Cluster = SeuratObj$Clusters)
    if(is.null(normal_groups)){
      tmpdat <- FetchData(SeuratObj, vars = c("Clusters", "EPCAM", "CDH1")) %>%
          mutate(Epi = EPCAM>0|CDH1>0) %>% group_by(Clusters) %>%
          summarize(Fraction = sum(Epi)/length(Epi)) %>%arrange(Fraction)
      normal_groups <- tmpdat$Clusters[tmpdat$Fraction<0.1]
      if(length(normal_groups)<2) normal_groups = tmpdat$Clusters[1:2]
    }
    if(is.null(gene_order_file)){
      gene_order_file <- system.file("extdata", "geneorder_hg20.txt.gz", package="scAnalyzer")
      gene_order_file <- read.table(gzfile(gene_order_file, "r"), header=FALSE,
                                    row.names=1, sep="\t", check.names=FALSE)
    }
    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = GetAssayData(SeuratObj, "counts"),
                                         gene_order_file = gene_order_file,
                                         annotations_file = ann,
                                         ref_group_names = normal_groups)
    infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=outdir,
                                  num_threads = nthreads,
                                  cluster_by_groups=TRUE, denoise=TRUE, HMM=TRUE)
  }
  if("numbat" %in% tolower(methods)){
    ## Add here
  }
  return(1)
}

