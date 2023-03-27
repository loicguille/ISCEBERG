#################################################################################################################################################################################################
############################################################################### READ FILE PAGE FUNCTIONS ########################################################################################
#################################################################################################################################################################################################

ReadCSV_files <- function(csvPath, metadataFiles = NULL, header = FALSE, separator, minGenebyCells, minCells){ # read CSV files
  tmp <- list()
  metadata <- list()
  if(header == "Yes"){
    head <- TRUE
  }else{
    head <- FALSE
  }
  if(separator == "tab"){
    sep <- "\t"
  }else{
    sep <- separator
  }
  
  if(is.null(metadataFiles)){
    for(i in 1:length(csvPath[,1])){
      tmp[[i]] <- read.csv(csvPath[[i,'datapath']],sep = sep, header = head)
      incProgress(1/(i+1), message = "Creating Seurat object")
      tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = minGenebyCells, min.features = minCells)
      tmp[[i]][["file"]] <- csvPath[[i,'name']]
    }
    seuratObj <- tmp[[1]]
    if(length(tmp) > 1){
      for(i in 2:length(csvPath[,1])){
        seuratObj <- merge(seuratObj,tmp[[i]])
      } 
    }
    seuratObj[["Project"]] <- "SeuratProject"
  }else{
    if(length(csvPath$datapath) != length(metadataFiles$datapath)){
      shinyalert("Oops", "The number of metadata files has to be the same as the number of input files", type = "error")
    }else{
      for(i in 1:length(csvPath[,1])){
        tmp[[i]] <- read.csv(csvPath[[i,'datapath']],sep = sep, header = headFile)
        metadata[[i]] <- read.csv(metadataFiles[[i,'datapath']])
        incProgress(1/(i+1), message = "Creating Seurat object")
        tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = minGenebyCells, min.features = minCells, meta.data = metadata[[i]])
        tmp[[i]][["file"]] <- csvPath[[i,'name']]
      }
      seuratObj <- tmp[[1]]
      if(length(tmp) > 1){
        for(i in 2:length(csvPath[,1])){
          seuratObj <- merge(seuratObj,tmp[[i]])
        } 
      }
      seuratObj[["Project"]] <- "SeuratProject"
    }
  }
  return(seuratObj)
}

readMatrix_files <- function(matrix_files, feature_files, cells_files, metadataFiles = NULL, minGenebyCells, minCells, FeatureName){### Read matrix files 
  tmp <- list()
  if(FeatureName == "EnsemblID"){
    col = 1
  }else{
    col = 2
  }
  if(length(matrix_files$datapath) != length(feature_files$datapath) || length(matrix_files$datapath) != length(cells_files$datapath)){
    shinyalert("Oops", "The number of matrix has to be the same as the number of cells files and feature files", type = "error")
    return(NULL)
  }else{
    for(i in 1:length(matrix_files[,1])){
      tmp[[i]] <- ReadMtx(mtx = matrix_files[[i,'datapath']], features = feature_files[[i,'datapath']], cells = cells_files[[i, 'datapath']], feature.column = col)
      incProgress(1/(i+1), message = "Creating Seurat object")
      tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = minGenebyCells, min.features = minCells)
      tmp[[i]][["file"]] <- matrix_files[[i,'name']]
    } 
    seuratObj <- tmp[[1]]
    if(length(tmp) > 1){
      for(i in 2:length(matrix_files[,1])){
        seuratObj <- merge(seuratObj,tmp[[i]])
      }   
    }
    seuratObj[["Project"]] <- "SeuratProject"
    return(seuratObj)
  }
}


readH5_files <- function(H5_files, metadataFiles = NULL, minGenebyCells, minCells){
  tmp <- list()
  metadata <- list()
  if(is.null(metadataFiles)){
    for(i in 1:length(H5_files[,1])){
      tmp[[i]] <- Read10X_h5(H5_files[[i,'datapath']])
      incProgress(1/(i+1), message = "Creating Seurat object")
      tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = minGenebyCells, min.features = minCells )
      tmp[[i]][["file"]] <- H5_files[[i,'name']]
    }
    seuratObj <- tmp[[1]]
    if(length(tmp) > 1){
      for(i in 2:length(H5_files[,1])){
        seuratObj <- merge(seuratObj,tmp[[i]])
      }
    }
    seuratObj[["Project"]] <- "SeuratProject"
    return(seuratObj)
  }else{
    if(length(H5_files$datapath) != length(metadataFiles$datapath)){
      shinyalert("Oops", "The number of metadata files has to be the same as the number of input files", type = "error")
    }else{
      for(i in 1:length(H5_files[,1])){
        tmp[[i]] <- Read10X_h5(H5_files[[i,'datapath']])
        metadata[[i]] <- read.csv(metadataFiles[[i,'datapath']])
        incProgress(1/(i+1), message = "Creating Seurat object")
        tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = minGenebyCells, min.features = minCells, meta.data = metadata[[i]])
        tmp[[i]][["file"]] <- H5_files[[i,'name']]
      }
      seuratObj <- tmp[[1]]
      if(length(tmp) > 1){
        for(i in 2:length(H5_files[,1])){
          seuratObj <- merge(seuratObj,tmp[[i]])
        } 
      }
      seuratObj[["Project"]] <- "SeuratProject"
      return(seuratObj)
    }
  }
}

Readtxt_files <- function(txt_files, metadataFiles = NULL, header = FALSE, separator, minGenebyCells, minCells){
  tmp <- list()
  metadata <- list()
  if(header=="Yes"){
    headFile <- TRUE
  }else{
    headFile <- FALSE
  }
  if(separator=="tab"){
    sep <- "\t"
  }else{
    sep <- separator
  }
  if(is.null(metadataFiles)){
    for(i in 1:length(txt_files[,1])){
      tmp[[i]] <- read.table(txt_files[[i,'datapath']], sep = sep, header = headFile)
      incProgress(1/(i+1), message = "Creating Seurat object")
      tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = minGenebyCells, min.features = minCells)
      tmp[[i]][["file"]] <- txt_files[[i,'name']]
    }
    seuratObj <- tmp[[1]]
    if(length(tmp) > 1){
      for(i in 2:length(txt_files[,1])){
        seuratObj <- merge(seuratObj,tmp[[i]])
      }   
    }
    seuratObj[["Project"]] <- "SeuratProject"
    return(seuratObj)
  }else{
    if(length(txt_files$datapath) != length(metadataFiles$datapath)){
      shinyalert("Oops", "The number of metadata files has to be the same as the number of input files", type = "error")
    }else{
      for(i in 1:length(txt_files[,1])){
        tmp[[i]] <- read.csv(txt_files[[i,'datapath']], sep = sep, header = headFile)
        metadata[[i]] <- read.csv(metadataFiles[[i,'datapath']])
        incProgress(1/(i+1), message = "Creating Seurat object")
        tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = minGenebyCells, min.features = minCells, meta.data = metadataFiles[[i]])
        tmp[[i]][["file"]] <- txt_files[[i,'name']]
      }
      seuratObj <- tmp[[1]]
      if(length(tmp) > 1){
        for(i in 2:length(txt_files[,1])){
          seuratObj <- merge(seuratObj,tmp[[i]])
        } 
      }
      seuratObj[["Project"]] <- "SeuratProject"
      return(seuratObj)
    }
  }
}
