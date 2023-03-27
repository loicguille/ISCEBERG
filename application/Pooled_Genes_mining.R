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

source("Utils.R")
poolGene <- function(obj, gene = NULL, typeOfNorm, coord, file = NULL, BoolColScale = FALSE, colorScale = NULL, size.points ){
  if(!is.null(file)){
    gene <- Readfile(obj, file)
  }
  if(typeOfNorm =="Sum" | typeOfNorm == "Mean"){
    if (typeOfNorm == "Sum"){
      pool_genes <- rowSums(FetchData(obj, vars =c(gene)))  
    }else if (typeOfNorm == "Mean"){
      pool_genes <- rowMeans(FetchData(obj, vars =c(gene)))
    }
    pool_genes <- cbind.data.frame(coord,counts = pool_genes)
    if(BoolColScale == TRUE){
      ggplot(pool_genes)+geom_point(aes(UMAP_1,UMAP_2, color = counts),shape=20, size=size.points)+ scale_color_gradient(low="lightgrey",high=colorScale) + theme_cowplot() 
    }else{
      ggplot(pool_genes)+geom_point(aes(UMAP_1,UMAP_2, color = counts),shape=20, size=size.points)+ scale_color_gradient(low="lightgrey",high="blue") + theme_cowplot()
    }
  }else{
    sbt <- AddModuleScore(obj, features =list(gene), name = "Gene_list")
    if(BoolColScale == TRUE){
      colorList <- c("lightgrey", colorScale)
      FeaturePlot(sbt, features = "Gene_list1", cols = colorList,pt.size = size.points)
    }else{
      FeaturePlot(sbt, features = "Gene_list1",pt.size = size.points)
    }
    
  }
  
}



VlnPlotPooled <- function(obj, gene = NULL, typeOfNorm, annotOrRes, cluster = NULL, annotation = NULL, file = NULL, color =NULL, color_list = NULL){
  if(!is.null(file)){
    gene <- Readfile(obj, file)
  }
  if(annotOrRes == "Resolution"){
    if(!is.null(color_list[[cluster]]$color)){
      color <- unlist(color_list[[cluster]]$color)
    }
  }else{
    if(!is.null(color_list[[annotation]]$color)){
      color <- unlist(color_list[[annotation]]$color)
    }
  }
  if(typeOfNorm =="Sum" | typeOfNorm == "Mean"){
    if(typeOfNorm == "Sum"){
      if(annotOrRes == "Resolution"){
        pool_genesVln <- cbind(obj[[cluster]],rowSums(FetchData(obj, vars = gene)))
      }else{
        pool_genesVln <- cbind(obj[[annotation]],rowSums(FetchData(obj, vars = gene)))
      }
    }else{
      if(annotOrRes == "Resolution"){
        pool_genesVln <- cbind(obj[[cluster]],rowMeans(FetchData(obj, vars = gene)))
      }else{
        pool_genesVln <- cbind(obj[[annotation]],rowMeans(FetchData(obj, vars = gene)))
      }
    }
    colnames(pool_genesVln) <- c("resolution","GeneList")
    if(is.null(color)){
      ggplot(pool_genesVln, aes(x=resolution,y=GeneList, fill = resolution))+geom_violin()+scale_fill_hue()+theme_cowplot()
    }else{
      ggplot(pool_genesVln, aes(x=resolution,y=GeneList, fill = resolution))+geom_violin()+scale_fill_hue()+theme_cowplot()+scale_fill_manual(values=color)
    }
  }else{
    sbt <- AddModuleScore(obj, features =list(gene), name = "Gene_list")
    if( annotOrRes == "Resolution"){
      VlnPlot(sbt, features = "Gene_list1", pt.size = 0, group.by = cluster, cols = color)
    }else{
      VlnPlot(sbt, features = "Gene_list1", pt.size = 0, group.by = annotation, cols = color)
    }
  }
}


DotplotPooled <- function(obj, gene = NULL, typeOfNorm, annotOrRes, cluster = NULL, annotation = NULL, file = NULL, BoolColScale = FALSE, colorScale = NULL){
  if(!is.null(file)){
    gene <- Readfile(obj, file)
  }
  if(typeOfNorm == "Sum"){
    if(annotOrRes == "Resolution"){
      Idents(obj) <- cluster
      test <- as.data.frame(cbind(rowSums(FetchData(obj, vars =c(gene))),obj[[cluster]]))
    }else{
      Idents(obj) <- annotation
      test <- as.data.frame(cbind(rowSums(FetchData(obj, vars =c(gene))),obj[[annotation]]))
    }
    table_percent <- CreateTableForPlot(test,obj, typeOfNorm = typeOfNorm)
    if(BoolColScale == TRUE){
      ggplot(table_percent , aes(x= GeneList, y  = cluster , size = percentExp, color = SumExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= colorScale)+theme_cowplot()
    }else{
      ggplot(table_percent , aes(x= GeneList, y  = cluster , size = percentExp, color = SumExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")+theme_cowplot()
    }
  }else if(typeOfNorm=="Mean"){
    if(annotOrRes == "Resolution"){
      Idents(obj) <- cluster
      test <- as.data.frame(cbind(rowMeans(FetchData(obj, vars =c(gene))),obj[[cluster]]))
    }else{
      Idents(obj) <- annotation
      test <- as.data.frame(cbind(rowMeans(FetchData(obj, vars =c(gene))),obj[[annotation]]))
    }
    table_percent <- CreateTableForPlot(test, obj,typeOfNorm = typeOfNorm )
    if(BoolColScale == TRUE){
      ggplot(table_percent , aes(x= GeneList, y  = cluster , size = percentExp, color = AvgExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= colorScale)+theme_cowplot()
    }else{
      ggplot(table_percent , aes(x= GeneList, y  = cluster , size = percentExp, color = AvgExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")+theme_cowplot()
    }
  }else{
    sbt <- AddModuleScore(obj, features =list(gene), name = "Gene_list")
    if(annotOrRes == "Resolution"){
      if(BoolColScale == TRUE){
        colorList <- c("lightgrey", colorScale)
        DotPlot(sbt, features = "Gene_list1", group.by = cluster, cols = colorList)
      }else{
        DotPlot(sbt, features = "Gene_list1", group.by = cluster)
      }
    }else{
      if(BoolColScale == TRUE){
        colorList <- c("lightgrey", colorScale)
        DotPlot(sbt, features = "Gene_list1", group.by = annotation, cols = colorList)
      }else{
        DotPlot(sbt, features = "Gene_list1", group.by = annotation)
      }
    }
  }
  
}

# Function used for creation of dotplot table (function DotplotPooled) Here is an example of what will contain percentExpressed at the end if you do a sum 
# cluster                                      logical number total_cells percentExp SumExp GeneList 
# <chr>                                        <lgl>    <int>       <int>      <dbl>  <dbl> <chr>    
# 1 GSE90047_BEC_Epcam+_dlk+_E17.5               TRUE        89          89      1     284.   Gene_list
# 2 GSE90047_Hep_Epcam+_dlk+_E17.5               TRUE        18          33      0.545   3.81 Gene_list
# 3 GSE90047_Hepatoblast_Epcam+_dlk+_E10.5_E13.5 TRUE       105         140      0.75   48.8  Gene_list
# 4 GSE90047_Hepatoblast_Epcam+_dlk+_E13.5_E15.5 TRUE        96         185      0.519  46.9  Gene_list
CreateTableForPlot <- function(table_test, obj, typeOfNorm){
  colnames(table_test)<- c("V1","V2")
  percentExpressed <- table_test %>% group_by(V2) %>%count(V1 > 0)
  colnames(percentExpressed) <- c("cluster", "logical","number")
  percentExpressed <- filter(percentExpressed, logical == "TRUE")
  totalCellsbyCluster <- table_test %>% group_by(V2) %>% tally()
  percentExpressed <- cbind(percentExpressed,total_cells = totalCellsbyCluster$n)
  percentExpressed$percentExp <- percentExpressed$number/percentExpressed$total_cells
  if(typeOfNorm == "Sum"){
    avgExpression <-table_test %>% group_by(V2) %>%summarise(avg =sum(V1))
    percentExpressed <- cbind(percentExpressed, SumExp =avgExpression$avg, GeneList = rep("Gene_list",length(percentExpressed$cluster)))
  }else{
    avgExpression <-table_test %>% group_by(V2) %>%summarise(avg =mean(V1))
    percentExpressed <- cbind(percentExpressed, AvgExp =avgExpression$avg, GeneList = rep("Gene_list",length(percentExpressed$cluster)))
  }
  return(percentExpressed)
}


