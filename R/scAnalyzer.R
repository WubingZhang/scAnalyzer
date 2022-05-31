#' Analysis pipeline for single cell RNA-seq data
#'
#' @docType methods
#' @name scAnalyzer
#' @rdname scAnalyzer
#'
#' @param obj A matrix or data.frame of count table or a Seurat object.
#' @param project A character, specifying the name of the project (also prefix of outputs).
#' @param analyses A character vector, specifying the type of analysis to perform,
#' available analysis include qc, doublet, norm, pca, clustering.
#' @param norm.method LogNormalize or SCTransform.
#' @param outdir The path to output figures.
#' @param resolution The resolution for clustering. 0.5 by default.
#' @param min.cells Include features detected in at least this many cells.
#' @param min.features Include cells with more than this many features are detected.
#' @param max.features Include cells with less than this many features are detected.
#' @param percent.mt Include cells with less than this fraction of mitochondrial gene expression.
#' @param subset.singlet Remove doublets from the data by subseting the singlets.
#' @param scale.factor Scale factor.
#' @param nVarfeatures The number of variable genes for downstream analysis.
#' @param nPCs The number of PCs.
#' @param PC.Variance At least this fraction of variance (80% by default) should be
#' explained by the Principal components, required when nPCs is not specified.
#'
#' @author Wubing Zhang
#'
#' @return Seurat object
#' @details This pipeline takes count table or seurat object as input, and
#' performs data quality control, normalization, doublets removal (DoubletFinder),
#' dimention reduction, cell clustering and identification of marker genes.
#'
#' @examples
#' library(scAnalyzer);
#' pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
#' res <- scAnalyzer(pbmc.data)
#' @import Seurat ggplot2 dplyr
#' @export
scAnalyzer <- function(obj, project = NULL,
                       analyses = c("qc", "doublet", "norm", "pca", "clustering"),
                       norm.method = c("LogNormalize", "SCTransform"),
                       outdir = "",
                       resolution = 0.5,
                       min.cells = 3,
                       min.features = 200,
                       max.features = 8000,
                       percent.mt = NULL,
                       subset.singlet = TRUE,
                       scale.factor = 1e4,
                       nVarfeatures = 2000,
                       nPCs = NULL,
                       PC.Variance = 0.8){
  requireNamespace("Seurat") || stop("Please install Seurat")
  requireNamespace("ggplot2") || stop("Please install ggplot2")
  requireNamespace("dplyr") || stop("Please install dplyr")
  #### Create Seurat Object ####
  if(class(obj)[1] != "Seurat"){
    message(Sys.time(), " Create Seurat Object ...")
    obj <- CreateSeuratObject(counts = obj, project = project,
                              min.cells = min.cells, min.features = min.features)
  }
  #### Data quality control ####
  if("qc" %in% tolower(analyses)){
    message(Sys.time(), " QC and normalization ...")
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^Mt-")
    if(dir.exists(outdir)){
      p_mt_1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      ggsave(plot = p_mt_1, paste0(outdir, "/VlnPlot_QC_1_", Sys.Date(), ".pdf"), width = 6, height = 4)
    }
    ## Filtering single cells with high mitochondrial gene expression
    if(!is.null(percent.mt)){
      obj <- subset(obj, subset = percent.mt < percent.mt)
    }else{
      obj <- subset(obj, subset = percent.mt < 50)
      if(sum(obj$percent.mt>5)/ncol(obj)<0.15){
        obj <- subset(obj, subset = percent.mt < 5)
      }else if(sum(obj$percent.mt>10)/ncol(obj)<0.15){
        if(mean(obj$percent.mt)+1.5*sd(obj$percent.mt)<10){
          obj <- subset(obj, subset = percent.mt < mean(percent.mt) + 1.5*sd(percent.mt))
        }else obj <- subset(obj, subset = percent.mt < 10)
      }else if(sum(obj$percent.mt>20)/ncol(obj)<0.15){
        if(mean(obj$percent.mt)+1.5*sd(obj$percent.mt)<20){
          obj <- subset(obj, subset = percent.mt < mean(percent.mt) + 1.5*sd(percent.mt))
        }else obj <- subset(obj, subset = percent.mt < 20)
      }else{
        obj <- subset(obj, subset = percent.mt < mean(percent.mt) + 1.5*sd(percent.mt))
        obj <- subset(obj, subset = percent.mt < 30)
      }
    }

    ## Filtering single cells with high total reads
    obj <- subset(obj, subset = nCount_RNA < mean(nCount_RNA) + 3*sd(nCount_RNA))
    obj <- subset(obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features)
    if(dir.exists(outdir)){
      p_mt_2 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      ggsave(plot = p_mt_2, paste0(outdir, "/VlnPlot_QC_2_", Sys.Date(), ".pdf"), width = 6, height = 4)
    }

    #### Remove doublets ####
    if("doublet" %in% tolower(analyses)){
      requireNamespace("DoubletFinder") || stop("Please install DoubletFinder")
      message(Sys.time(), " Predicting doublets using DoubletFinder ...")

      # Identify cell clusters
      if(tolower(norm.method[1])=="sctransform"){
        obj <- SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
      }else{
        obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = FALSE)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nVarfeatures)
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
      if(subset.singlet) obj <- subset(obj, subset = DF.class=="Singlet")
    }
  }

  #### Data normalization ####
  if("norm" %in% tolower(analyses)){
    message(Sys.time(), " Data normalization ...")
    if(tolower(norm.method[1])=="sctransform"){
      obj <- SCTransform(obj, method = "glmGamPoi", return.only.var.genes = FALSE, verbose = FALSE)
    }else{
      obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = scale.factor, verbose = FALSE)
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nVarfeatures)
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
    if(is.null(nPCs)){
      # Selecting PCs capturing more than 80% variance
      vars = Stdev(obj, reduction = "pca")^2
      PCs = seq(1, which(cumsum(vars / sum(vars)) > PC.Variance)[1])
    }else{
      PCs = seq(1, nPCs)
    }
    message("# of PCs for downstream analysis: ", max(PCs))
    obj <- FindNeighbors(object = obj, dims = PCs, verbose = FALSE)
    ## Using Leiden algorithm to detect the cell communities (leidenalg python)
    obj <- FindClusters(object = obj, resolution = resolution, algorithm = 1, verbose = FALSE)
    obj <- RunUMAP(obj, dims = PCs, verbose = FALSE)
    if(dir.exists(outdir)){
      p_umap <- DimPlot(obj, reduction = "umap", label=TRUE, label.size = 8)
      ggsave(plot = p_umap, paste0(outdir, "/UmapPlot_", Sys.Date(), ".pdf"), width = 6, height = 5)
    }
  }

  if(dir.exists(outdir)){
    saveRDS(obj, paste0(outdir, "/", project, "_scRNA_obj.rds"))
  }
  return(obj)
}
