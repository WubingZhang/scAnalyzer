#' Analysis pipeline for single cell RNA-seq data
#'
#' @docType methods
#' @name scAnalyzer
#' @rdname scAnalyzer
#'
#' @param obj A matrix or data.frame of count table or a seurat object.
#' @param project A character, specifying the name of the project, which will be used in output file names.
#' @param analyses A character vector, specifying the type of analysis to perform, avaible analysis include qc,
#' doublet, norm, pca, clustering, findmarker.
#' @param norm.method LogNormalize or SCTransform.
#' @param resolution The resolution for clustering. 0.1 by default.
#' @param outdir The path to output figures.
#'
#' @author Wubing Zhang
#'
#' @return A list a seurat object and markers for all clusters
#' @details This pipeline takes count table as input, and performs data quality control, normalization, doublets removal (scDblFinder),
#' dimention reduction, cell clustering and identification of marker genes.
#'
#' @examples
#' library(scAnalyzer);
#' pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' res <- scAnalyzer(pbmc.data)
#' @import Seurat DoubletFinder ggplot2 dplyr
#' @export

scAnalyzer <- function(obj, project = NULL,
                       analyses = c("qc", "doublet", "norm", "pca", "clustering"),
                       norm.method = c("LogNormalize", "SCTransform"),
                       resolution = 0.5,
                       outdir = ""){
  requireNamespace("Seurat") || stop("Please install Seurat")
  requireNamespace("ggplot2") || stop("Please install ggplot2")
  requireNamespace("dplyr") || stop("Please install dplyr")
  #### Create Seurat Object ####
  if(class(obj) != "Seurat"){
    message(Sys.time(), " Create Seurat Object ...")
    obj <- CreateSeuratObject(counts = obj, project = project, min.cells = 3, min.features = 200)
  }
  #### Data quality control ####
  if("qc" %in% tolower(analyses)){
    message(Sys.time(), " QC and normalization ...")
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    if(dir.exists(outdir)){
      p_mt_1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      ggsave(plot = p_mt_1, paste0(outdir, "/VlnPlot_QC_1_", Sys.Date(), ".pdf"), width = 6, height = 4)
    }
    ## Filtering single cells with high mitochondrial gene expression
    obj <- subset(obj, subset = percent.mt < 50)
    if(sum(obj$percent.mt>5)/ncol(obj)<0.2){
      obj <- subset(obj, subset = percent.mt < 5)
    }else if(sum(obj$percent.mt>10)/ncol(obj)<0.2){
      if(mean(obj$percent.mt)+1.5*sd(obj$percent.mt)<10){
        obj <- subset(obj, subset = percent.mt < mean(percent.mt) + 1.5*sd(percent.mt))
      }else obj <- subset(obj, subset = percent.mt < 10)
    }else if(sum(obj$percent.mt>20)/ncol(obj)<0.2){
      if(mean(obj$percent.mt)+1.5*sd(obj$percent.mt)<20){
        obj <- subset(obj, subset = percent.mt < mean(percent.mt) + 1.5*sd(percent.mt))
      }else obj <- subset(obj, subset = percent.mt < 20)
    }else{
      obj <- subset(obj, subset = percent.mt < mean(percent.mt) + 1.5*sd(percent.mt))
      obj <- subset(obj, subset = percent.mt < 30)
    }
    ## Filtering single cells with high total reads
    obj <- subset(obj, subset = nCount_RNA < mean(nCount_RNA) + 3*sd(nCount_RNA))
    obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000)
    if(dir.exists(outdir)){
      p_mt_2 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      ggsave(plot = p_mt_2, paste0(outdir, "/VlnPlot_QC_2_", Sys.Date(), ".pdf"), width = 6, height = 4)
    }

    #### Remove doublets ####
    if("doublet" %in% tolower(analyses)){
      requireNamespace("DoubletFinder") || stop("Please install DoubletFinder")
      message(Sys.time(), " Remove doublets ...")

      # Identify cell clusters
      if(tolower(norm.method[1])=="sctransform"){
        obj <- SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
      }else{
        obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
        obj <- ScaleData(obj)
      }
      obj <- RunPCA(object = obj, features = VariableFeatures(object = obj), verbose = FALSE)
      obj <- FindNeighbors(object = obj, dims = 1:10, verbose = FALSE)
      obj <- FindClusters(object = obj, resolution = 0.1, algorithm = 1, verbose = FALSE)

      ## Removing doublets using doubletFinder
      prints <- capture.output(suppressMessages(
        sweep.res.list <- paramSweep_v3(obj, PCs = 1:10, sct = FALSE, num.cores = 4)))
      prints <- capture.output(suppressMessages(
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)))
      bcmvn <- find.pK(sweep.stats)
      ## Homotypic Doublet Proportion Estimate
      annotations <- obj@meta.data$seurat_clusters
      homotypic.prop <- modelHomotypic(annotations)
      ## Estimating doublet formation rate
      # https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php
      nExp_poi <- round(nrow(obj@meta.data)*nrow(obj@meta.data)*7.75e-6)
      nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

      ## Run DoubletFinder with varying classification stringencies -----
      prints <- capture.output(suppressMessages(
        obj <- doubletFinder_v3(obj, PCs = 1:10, pK = 0.09, sct = FALSE, nExp = nExp_poi.adj)))
      colnames(obj@meta.data)[grepl("DF.classification", colnames(obj@meta.data))] <- "DF.class"
      colnames(obj@meta.data)[grepl("^pANN_", colnames(obj@meta.data))] <- "DoubletFinder"
      if(dir.exists(outdir)){
        p.obj <- DimPlot(obj, group.by = "DF.class")
        ggsave(plot = p.obj, paste0(outdir, "/DimPlot_Doublet.class_", Sys.Date(), ".pdf"), width = 4.5, height = 4)
        p_dbl = FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "DoubletFinder")
        ggsave(plot = p_dbl, paste0(outdir, "/FeatureScatter_nCount_Doublet_", Sys.Date(), ".pdf"), width = 4.5, height = 4)
        p_nfeat = FeaturePlot(obj, features = c("nCount_RNA", "DoubletFinder"))
        ggsave(plot = p_nfeat, paste0(outdir, "/FeaturePlot_nCount_Doublet_", Sys.Date(), ".pdf"), width = 7, height = 4)
      }
      obj <- subset(obj, subset = DF.class=="Singlet")
    }
  }

  #### Data normalization ####
  if("norm" %in% tolower(analyses)){
    message(Sys.time(), " Data normalization ...")
    if(tolower(norm.method[1])=="sctransform"){
      obj <- SCTransform(obj, method = "glmGamPoi", return.only.var.genes = FALSE, verbose = FALSE)
    }else{
      obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
      obj <- ScaleData(obj, features = rownames(obj))
    }
    if(dir.exists(outdir)){
      p_vargene <- LabelPoints(plot = VariableFeaturePlot(obj), points = VariableFeatures(obj)[1:10], repel = TRUE)
      ggsave(plot = p_vargene, paste0(outdir, "/VarFeaturePlot_", Sys.Date(), ".pdf"), width = 8, height = 4)
    }
  }

  #### Dimention reduction ####
  if("pca" %in% tolower(analyses)){
    message(Sys.time(), " Dimention reduction ...")
    obj <- RunPCA(object = obj, features = VariableFeatures(object = obj), verbose = FALSE)
    if(dir.exists(outdir)){
      p_pc1 <- VizDimLoadings(obj, dims = 1:2, reduction = "pca")
      ggsave(plot = p_pc1, paste0(outdir, "/VizDimLoadingsPlot_PC1_PC2_", Sys.Date(), ".pdf"), width = 6.5, height = 4.5)
      p_pc2 <- DimHeatmap(obj, dims = 1:9, cells = 500, balanced = FALSE, fast = FALSE)
      ggsave(plot = p_pc2, paste0(outdir, "/DimHeatmap_PC1-9_", Sys.Date(), ".pdf"), width = 11, height = 12)
      p_pc3 <- ElbowPlot(obj)
      ggsave(plot = p_pc3, paste0(outdir, "/ElbowPlot_", Sys.Date(), ".pdf"), width = 4.5, height = 3.5)
      # p_mito = FeaturePlot(obj, features = c("nFeature_RNA", "percent.mt"))
      # ggsave(plot = p_mito, paste0(outdir, "/FeaturePlot_QC_", Sys.Date(), ".pdf"), width = 7, height = 4)
    }
  }

  #### Clustering analysis ####
  if("clustering" %in% tolower(analyses)){
    message(Sys.time(), " Clustering analysis ...")
    # Selecting PCs capturing more than 80% variance
    vars = Stdev(obj, reduction = "pca")^2
    PCs = seq(1, which(cumsum(vars / sum(vars)) > 0.8)[1])
    message("# of PCs for downstream analysis: ", max(PCs))
    obj <- FindNeighbors(object = obj, dims = PCs, verbose = FALSE)
    ## Using Leiden algorithm to detect the cell communities (leidenalg python)
    obj <- FindClusters(object = obj, resolution = resolution, algorithm = 1, verbose = FALSE)
    obj <- RunUMAP(obj, dims = PCs, verbose = FALSE)
    if(dir.exists(outdir)){
      # p_nfeat = FeaturePlot(obj, features = c("nFeature_RNA", "percent.mt"))
      # ggsave(plot = p_nfeat, paste0(outdir, "/FeaturePlot_nFeature_mtPercent_", Sys.Date(), ".pdf"), width = 12, height = 3.8)
      p_umap <- DimPlot(obj, reduction = "umap", label=TRUE, label.size = 8)
      ggsave(plot = p_umap, paste0(outdir, "/UmapPlot_", Sys.Date(), ".pdf"), width = 6, height = 5)
    }
  }

  #### Identifying DE genes ####
  markers <- NULL
  if("findmarker" %in% tolower(analyses)){
    message(Sys.time(), " Differential expression analysis ...")
    markers<- presto::wilcoxauc(obj, 'seurat_clusters', assay = 'data')
    top5 <- top_markers(markers, n = 5, auc_min = .6, pct_in_min = 60)[,-1] %>% unlist() %>% unique()
    top10 <- top_markers(markers, n = 10, auc_min = .6, pct_in_min = 60)[,-1] %>% unlist() %>% unique()
    if(dir.exists(outdir)){
      saveRDS(markers, paste0(outdir, "/Seurat_", project, "_AllMarkers_", Sys.Date(), ".rds"))
      p1 <- DoHeatmap(obj, features = top5) + NoLegend()
      ggsave(plot = p1, paste0(outdir, "/Heatmap_Markers_top5_", Sys.Date(), ".pdf"),
             width = 8, height = 4+length(top5)/10)
      p2 <- DoHeatmap(obj, features = top10) + NoLegend()
      ggsave(plot = p1, paste0(outdir, "/Heatmap_Markers_top10_", Sys.Date(), ".pdf"),
             width = 8, height = 4+length(top10)/10)
    }
  }

  if(dir.exists(outdir)){
    saveRDS(obj, paste0(outdir, "/SeuratObj_", project, "_", Sys.Date(), ".rds"))
  }
  return(obj)
}
