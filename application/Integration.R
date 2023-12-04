#####################################################################################################
######################################### Integration ###############################################
#####################################################################################################
library(harmony)
library(STACAS)
library(scIntegrationMetrics)
library(tidyr)

harmony_integration <- function(object,normalization, varReg, minRes, maxRes, stepRes){
  object <- RunPCA(object,reduction.name = "pca_harmony", npcs = 30)
  object <- RunHarmony(object,reduction.use = "pca_harmony", group.by.vars = varReg, max.iter.harmony = 10)
  object <- RunUMAP(object, reduction = "harmony",dims = 1:30, reduction.name = "umap_harmony")
  object <- FindNeighbors(object, reduction = "harmony" ,dims = 1:30, graph.name = "Harmony_snn")
  for(i in seq(minRes, maxRes, stepRes)){
    object <- FindClusters(object, res = i, graph.name = "Harmony_snn")
  }
  return(object)
}

Seurat_integration <- function(object, normalization, minRes, maxRes, StepRes, reductionUsed, minGenes, maxGenes, maxCounts, maxPercent){
  object <- subset(object, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent & nCount_RNA < maxCounts)
  withProgress(message = "Running normalization", value = 1/4,{
    objectList <- SplitObject(object, split.by = "file")
    objectList <- lapply(X = objectList, FUN = function(x) {
      if(normalization == "SCTransform"){
        x <- SCTransform(x)
      }else{
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      }
    })
    all.genes_RNA <- row.names(objectList[[1]])
    for(i in 2 : length(objectList)){
      all.genes_RNA <- intersect(all.genes_RNA, rownames(objectList[[i]]))
    }
  })
  withProgress(message = "Select integration features", value = 2/4,{
    features <- SelectIntegrationFeatures(object.list = objectList)
  })
  withProgress(message = "Find anchors and integrate data", value = 3/4,{
    if(normalization == "SCTransform"){
      objectList <- PrepSCTIntegration(object.list = objectList, anchor.features = features)
      immune.anchors <- FindIntegrationAnchors(object.list = objectList, anchor.features = features, normalization.method = "SCT",reduction = reductionUsed)
      integration <- IntegrateData(anchorset = immune.anchors,normalization.method = "SCT",features.to.integrate = all.genes_RNA)
    }else{
      immune.anchors <- FindIntegrationAnchors(object.list = objectList, anchor.features = features,normalization.method = "LogNormalize",reduction = reductionUsed)
      integration <- IntegrateData(anchorset = immune.anchors,normalization.method = "LogNormalize",features.to.integrate = all.genes_RNA)
      integration <- ScaleData(integration)
    }
  })
  withProgress(message = "Calculate PCA and UMAP coordinate and clusterize dataset", value = 3.90/4,{
    integration <- RunPCA(integration,reduction.name = "pca_seurat", npcs = 30)
    integration <- RunUMAP(integration, reduction = "pca_seurat",dims = 1:30, reduction.name = "umap_seurat")
    integration <- FindNeighbors(integration, reduction = "pca_seurat", dims = 1:30, graph.name = "Seurat_snn")
    for(i in seq(minRes, maxRes, StepRes)){
      integration <- FindClusters(integration, res = i, graph.name = "Seurat_snn")
    }
  })
  # DefaultAssay(integration) <- "RNA"
  # if(normalization == "SCTransform"){
  #   integration <- SCTransform(integration)
  # }else{
  #   integration <- NormalizeData(integration)
  #   integration <- FindVariableFeatures(integration, selection.method = "vst", nfeatures = 2000)
  #   integration <- ScaleData(integration)
  # }
  # integration <- RunPCA(integration, npcs = 30, reduction.name ="pca")
  # integration <- RunUMAP(integration, reduction = "pca",dims = 1:30, reduction.name = "umap")
  # integration <- FindNeighbors(integration, reduction = "pca", dims = 1:30)
  # for(i in seq(minRes, maxRes, stepRes)){
  #   integration <- FindClusters(integration, res = i)
  # }
  return(integration)
}


STACAS_integration <- function(object, minGenes, maxGenes, maxCounts, maxPercent, minRes, maxRes, stepRes){
  object <- subset(object, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent & nCount_RNA < maxCounts)
  objectList <- SplitObject(object, split.by = "file")
  objectList <- lapply(X = objectList, FUN = function(x){
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  all.genes_RNA <- row.names(objectList[[1]])
  for(i in 2 : length(objectList)){
    all.genes_RNA <- intersect(all.genes_RNA, rownames(objectList[[i]]))
  }
 
  features <- SelectIntegrationFeatures(object.list = objectList, nfeatures = 2000)
  Anchors_found <- FindAnchors.STACAS(object.list = objectList, dims=1:30, anchor.features=features, alpha = 0.8, anchor.coverage = 0.5, verbose = T)
  Integrated_object <- IntegrateData.STACAS(Anchors_found, dims=1:30, features.to.integrate=all.genes_RNA, semisupervised = FALSE)
  Integrated_object <- ScaleData(Integrated_object)
  Integrated_object <- RunPCA(Integrated_object, verbose = FALSE, npcs = 30, reduction.name = "pca_stacas")
  Integrated_object <- RunUMAP(Integrated_object, reduction = "pca_stacas", dims = 1:30, reduction.name = "umap_stacas")
  Integrated_object <- FindNeighbors(Integrated_object, reduction = "pca_stacas", graph.name = "STACAS_snn")
  for(i in seq(minRes, maxRes,stepRes)){
    Integrated_object <- FindClusters(Integrated_object, res = i, graph.name = "STACAS_snn")
  }
  # DefaultAssay(Integrated_object) <- "RNA"
  # Integrated_object <- NormalizeData(Integrated_object)
  # Integrated_object <- FindVariableFeatures(Integrated_object, selection.method = "vst", nfeatures = 2000)
  # Integrated_object <- ScaleData(Integrated_object)
  # Integrated_object <- RunPCA(Integrated_object, npcs = 30, reduction.name ="pca")
  # Integrated_object <- RunUMAP(Integrated_object, reduction = "pca",dims = 1:30, reduction.name = "umap")
  # Integrated_object <- FindNeighbors(Integrated_object, reduction = "pca", dims = 1:30)
  # for(i in seq(minRes, maxRes, stepRes)){
  #   Integrated_object <- FindClusters(Integrated_object, res = i)
  # }
  return(Integrated_object)
}

Integration_according_to_interest <- function(object, method, normalization, minGenes, maxGenes, maxCounts, maxPercent, minRes, maxRes, StepRes, reductionUsed = NULL, varReg = NULL){
  integration <- subset(object, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent & nCount_RNA < maxCounts)
  objectList <- SplitObject(integration, split.by = "file")
  objectList <- lapply(X = objectList, FUN = function(x) {
    if(normalization == "SCTransform"){
      x <- SCTransform(x)
    }else{
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    }
  })
  all.genes_RNA <- row.names(objectList[[1]])
  for(i in 2 : length(objectList)){
    all.genes_RNA <- intersect(all.genes_RNA, rownames(objectList[[i]]))
  }
  if("Seurat" %in% method){
    features <- SelectIntegrationFeatures(object.list = objectList, nfeatures = 2000)
    if(normalization == "SCTransform"){
      objectList <- PrepSCTIntegration(object.list = objectList, anchor.features = features)
      if(reductionUsed == "rpca"){
        objectList <- lapply(X = objectList, FUN = function(x) {
          x <- RunPCA(x)
        })
      }
      immune.anchors <- FindIntegrationAnchors(object.list = objectList, anchor.features = features, normalization.method = "SCT",reduction = reductionUsed)
      integration <- IntegrateData(anchorset = immune.anchors,normalization.method = "SCT",features.to.integrate = all.genes_RNA)
    }else{
      if(reductionUsed == "rpca"){
        objectList <- lapply(X = objectList, FUN = function(x) {
          x <- ScaleData(x)
          x <- RunPCA(x)
        })
      }
      immune.anchors <- FindIntegrationAnchors(object.list = objectList, anchor.features = features,normalization.method = "LogNormalize",reduction = reductionUsed)
      integration <- IntegrateData(anchorset = immune.anchors,normalization.method = "LogNormalize",features.to.integrate = all.genes_RNA)
      integration <- ScaleData(integration)
    }
    integration <- RunPCA(integration,reduction.name = "pca_seurat", npcs = 30)
    integration <- RunUMAP(integration, reduction = "pca_seurat",dims = 1:30, reduction.name = "umap_seurat")
    integration <- FindNeighbors(integration, reduction = "pca_seurat", dims = 1:30, graph.name = "Seurat_snn")
    for(i in seq(minRes, maxRes, StepRes)){
      integration <- FindClusters(integration, res = i, graph.name = "Seurat_snn")
    }
  }
  if("STACAS" %in% method){
    features <- SelectIntegrationFeatures(object.list = objectList, nfeatures = 2000)
    Anchors_found <- FindAnchors.STACAS(object.list = objectList, dims=1:30, anchor.features=features, alpha = 0.8, anchor.coverage = 0.5, verbose = T)
    integration <- IntegrateData.STACAS(Anchors_found, dims=1:30, features.to.integrate=all.genes_RNA, semisupervised = FALSE, new.assay.name = "integrated_stacas")
    integration <- ScaleData(integration)
    integration <- RunPCA(integration, verbose = FALSE, npcs = 30, reduction.name = "pca_stacas")
    integration <- RunUMAP(integration, reduction = "pca_stacas", dims = 1:30, reduction.name = "umap_stacas")
    integration <- FindNeighbors(integration, reduction = "pca_stacas", graph.name = "STACAS_snn")
    for(i in seq(minRes, maxRes,StepRes)){
      integration <- FindClusters(integration, res = i, graph.name = "STACAS_snn")
    }
  }
  DefaultAssay(integration) <- "RNA"
  if(normalization == "SCTransform"){
    integration <- SCTransform(integration)
  }else{
    integration <- NormalizeData(integration)
    integration <- FindVariableFeatures(integration)
    integration <- ScaleData(integration)
  }
  integration <- RunPCA(integration)
  integration <- FindNeighbors(integration,dims = 1:30) 
  integration <- RunUMAP(integration, dims = 1:30) 
  for(i in seq(minRes,maxRes, StepRes)){
    integration <- FindClusters(integration, res = i)
  }
  if("harmony"%in% method){
    integration <- RunPCA(integration,reduction.name = "pca_harmony", npcs = 30)
    integration <- RunHarmony(integration,reduction.use = "pca_harmony", group.by.vars = varReg, max.iter.harmony = 10)
    integration <- RunUMAP(integration, reduction = "harmony",dims = 1:30, reduction.name = "umap_harmony")
    integration <- FindNeighbors(integration, reduction = "harmony" ,dims = 1:30, graph.name = "Harmony_snn")
    for(i in seq(minRes, maxRes, StepRes)){
      integration <- FindClusters(integration, res = i, graph.name = "Harmony_snn")
    }
  }
  return(integration)
}

CreatePlotIntegration <- function(object,batch,method){
  # batch value contain information about batch in the dataset
  LISI <- list()
  lisi_perplexity <- 30
  
  ## For merged datasets
  method_name <- "Merged"
  LISI[[method_name]] <- compute_lisi(X = object@reductions[["umap"]]@cell.embeddings, meta_data = object@meta.data, label_colnames=batch, perplexity = lisi_perplexity)
  names(LISI$Merged) <- "batch"
  lisi <- data.frame("Merged" = LISI$Merged$batch)
  
  ## For integrated dataset
  if("harmony" %in% method){
    method_name <- "harmony"
    LISI[[method_name]]<- compute_lisi(X = object@reductions[["umap_harmony"]]@cell.embeddings, meta_data = object@meta.data, label_colnames=batch, perplexity = lisi_perplexity)
    names(LISI$harmony) <- "batch"
    lisi <- cbind(lisi,data.frame("harmony" = LISI$harmony$batch))
  }
  if("Seurat" %in% method){
    method_name <- "Seurat"
    LISI[[method_name]]<- compute_lisi(X = object@reductions[["umap_seurat"]]@cell.embeddings, meta_data = object@meta.data, label_colnames=batch, perplexity = lisi_perplexity)
    names(LISI$Seurat) <- "batch"
    lisi <- cbind(lisi,data.frame("Seurat" = LISI$Seurat$batch))
  }
  if("STACAS" %in% method){
    method_name <- "STACAS"
    LISI[[method_name]]<- compute_lisi(X = object@reductions[["umap_stacas"]]@cell.embeddings, meta_data = object@meta.data, label_colnames=batch, perplexity = lisi_perplexity)
    names(LISI$STACAS) <- "batch"
    lisi <- cbind(lisi,data.frame("STACAS" = LISI$STACAS$batch))
  }
  lisi <- gather(lisi, Methods, LisiValue)
  ggplot(lisi)+geom_boxplot(aes(x=Methods, y = LisiValue, fill = Methods ))+ggtitle("Batch LISI for different integration method")+theme_cowplot()
}

CreatePlotSilhouette <- function(object, bioconservation){
  # bioconservation value contain information about cell type in the dataset
  Silhouette <- list()
  lisi_perplexity <- 30
  
  method <- "Integrated"
  Silhouette[[method]] <- compute_silhouette(X = object@reductions[["umap"]]@cell.embeddings, meta_data = object@meta.data, label_colnames = bioconservation)
  names(Silhouette$Integrated) <- "bioconservation"
  conservation <- data.frame("Integrated" = Silhouette$Integrated$bioconservation)
  conservation <- gather(conservation, Methods, LisiValue)
  ggplot(conservation)+geom_boxplot(aes(x=Methods, y = LisiValue, fill = Methods ))+ggtitle("Bioconservation for different integration method")
  
}