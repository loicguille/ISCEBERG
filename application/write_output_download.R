write_color_file<- function(color_list, name_file){
  colors <- sapply(color_list, function(x) x$color)
  labels <- sapply(color_list, function(x) x$label)
  
  max_length <- max(sapply(colors, length))
  
  
  for (i in seq_along(colors)) {
    if (length(colors[[i]]) < max_length) {
      colors[[i]] <- c(colors[[i]], rep(NA, max_length - length(colors[[i]])))
      labels[[i]] <- c(labels[[i]], rep(NA, max_length - length(labels[[i]])))
    }
  }
  df <- data.frame(colors)
  names(df) <- paste0(names(color_list), ".color")
  dt <- data.frame(labels)
  names(dt) <- paste0(names(color_list), ".annot")
  dataFrameColLabel <- cbind(dt,df)
  dataFrameColLabel<- dataFrameColLabel[, order(names(dataFrameColLabel))]
  write.csv(dataFrameColLabel, name_file, row.names = FALSE,)
}

write_rmarkdown_report_preprocessing <- function(name_file, normalization, minGenes, maxGenes, maxCount, maxPercent,resMax, resMin, resStep){
  tempReport <- file.path(tempdir(),"logFile.Rmd")
  file.copy("logFile.Rmd",tempReport ,overwrite = TRUE)
  command <- list()
  line <- paste0("SeuratObjsubset <- subset(seuratObj, subset = nFeature_RNA > ",minGenes," & nFeature_RNA < ",maxGenes," & nCount_RNA < ",maxCount," & percent.mt < ",maxPercent,")")
  command <- append(command,line)
  
  if(normalization == "SCTransform"){
    line <- paste0("SeuratObjsubset <- SCTransform(SeuratObjsubset)")
    command <- append(command,line)
  }else{
    line <- paste0("SeuratObjsubset <- NormalizeData(SeuratObjsubset)")
    command <- append(command,line)
    line <- paste0("SeuratObjsubset <- FindVariableFeatures(SeuratObjsubset)")
    command <- append(command,line)
    line <- paste0("SeuratObjsubset <- ScaleData(SeuratObjsubset)")
    command <- append(command,line)
  }
  command <- append(command,"SeuratObjsubset <- RunPCA(SeuratObjsubset)")
  command <- append(command,"SeuratObjsubset <- FindNeighbors(SeuratObjsubset, dims  = 1:30)")
  command <- append(command,"SeuratObjsubset <- RunUMAP(SeuratObjsubset, dims = 1:30)")
  for(i in seq(resMin,resMax, resStep)){
    command <- append(command,paste0("SeuratObjsubset <- FindClusters(SeuratObjsubset, res =", i,")"))
  }
  params <- list(use = command)
  rmarkdown::render(tempReport,output_file = name_file,params = params, envir = new.env(parent = globalenv()))
}


write_rmarkdown_report_subclustering <- function(name_file, ClustOrAnnot,cluster, subsetCluster, Annot, subsetAnnot){
  tempReport <- file.path(tempdir(),"logFile.Rmd")
  file.copy("logFile.Rmd",tempReport ,overwrite = TRUE)
  command <- list()
  if(ClustOrAnnot == "Cluster"){
    line <- paste0("subsetSeuratObj <- subset(seuratObj, subset =",cluster," == ",list(subsetCluster),")")
  }else{
    line <- paste0("subsetSeuratObj <- subset(seuratObj, subset = ",Annot," == ",list(subsetAnnot),")")
  }
  command <- append(command,line)
  command <- append(command,"SeuratObjsubset <- FindVariableFeatures(SeuratObjsubset)")
  command <- append(command,"SeuratObjsubset <- RunPCA(SeuratObjsubset)")
  command <- append(command,"SeuratObjsubset <- FindNeighbors(SeuratObjsubset)")
  command <- append(command,"SeuratObjsubset <- RunUMAP(SeuratObjsubset, dims = 1:30)")
  for(i in seq(0,1,0.1)){
    command <- append(command,paste0("SeuratObjsubset <- FindClusters(SeuratObjsubset, res =", i,")"))
  }
  params <- list(use = command)
  rmarkdown::render(tempReport,output_file = name_file,params = params, envir = new.env(parent = globalenv()))
}
