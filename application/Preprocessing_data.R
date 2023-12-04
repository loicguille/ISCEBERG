library(Seurat)
Preprocessing_seuratObject <- function(object, normalization, minGenes, maxGenes, maxCount, maxPercent,resMax, resMin, resStep){
  object <- subset(object, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent & nCount_RNA < maxCount)
  DefaultAssay(object) <- "RNA"
  withProgress(message = "Normalization",value =1/3,{
    if(normalization == "SCTransform"){
      object <- SCTransform(object)
    }else{
      object <- NormalizeData(object)
      object <- FindVariableFeatures(object)
      object <- ScaleData(object)
    }
  })
  withProgress(message ="Calculate PCA and UMAP coordinate", value = 2/3,{
    object <- RunPCA(object)
    object <- FindNeighbors(object,dims = 1:30) 
    object <- RunUMAP(object, dims = 1:30) 
  })
  withProgress(message ="Clusterization", value = 2.90/3,{
    for(i in seq(resMin,resMax, resStep)){
      object <- FindClusters(object, res = i)
    }
  })
  return(object)
}
