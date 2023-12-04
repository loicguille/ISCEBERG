### function for annotation page
library(Seurat)
library(Matrix)
library(bslib)
library(ggplot2)
library(DT)
library(cowplot)
library(RColorBrewer)
library(plotly)
library(grid)
library(dplyr)
library(gridExtra)

#obj = seurat object, nameVar = name of existing annotation, annotation = text to put in annotation, resolution = chosen resolution for annotation
#clusterToAnnotate = which cluster(s) will be annotated

Annotation_existent <- function(obj, nameVar, annotation, resolution, clusterToAnnotate){
  f.clusters <- obj[[nameVar]]
  f.clusters[obj[[resolution]] == clusterToAnnotate[1]] <- annotation
  if(length(clusterToAnnotate) > 1){
    for(i in 2:length(clusterToAnnotate)){
      f.clusters[obj[[resolution]] == clusterToAnnotate[i]] <- annotation
    }
  }
  return(f.clusters)
}


Annotation_unexistent <- function(obj, annotation, resolution, clusterToAnnotate){
  f.clusters <- rep("unknown", dim(obj)[2])
  names(f.clusters) <- colnames(obj)
  f.clusters[obj[[resolution]] == clusterToAnnotate[1]] <- annotation
  if(length(clusterToAnnotate) > 1){
    for(i in 2:length(clusterToAnnotate)){
      f.clusters[obj[[resolution]] == clusterToAnnotate[i]] <- annotation
    }
  }
  return(f.clusters)
}
